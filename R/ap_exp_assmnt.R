getRasterValueInRange_LargeScale <- function(
    dataframe,
    input_folder_name,
    output_folder_name,
    pollutant = "PM",
    use_parallel = TRUE,
    n_cores = NULL,
    batch_size = 50000,
    temp_dir = NULL
) {
    required_packages <- c("data.table", "terra", "lubridate")
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    if (length(missing_packages) > 0) {
        stop(sprintf("Required packages not installed: %s\nInstall with: install.packages(c('%s'))",
            paste(missing_packages, collapse = ", "),
            paste(missing_packages, collapse = "', '")))
    }
    library(data.table); library(terra); library(lubridate)

    if (use_parallel) {
        if (!requireNamespace("parallel", quietly = TRUE)) {
            warning("Package 'parallel' not available. Running in serial mode.")
            use_parallel <- FALSE
        } else {
            library(parallel)
            if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
            cat(sprintf("Parallel processing enabled with %d cores\n", n_cores))
        }
    }

    start_time <- Sys.time()
    required_cols <- c("lat", "long", "start_date", "end_date")
    missing_cols <- required_cols[!required_cols %in% names(dataframe)]
    if (length(missing_cols) > 0) stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")))
    if (any(is.na(dataframe$lat)) || any(is.na(dataframe$long))) stop("Missing lat/long in dataframe")

    # Ensure dates are correct
    dataframe$start_date <- as.Date(dataframe$start_date)
    dataframe$end_date   <- as.Date(dataframe$end_date)
    dt <- as.data.table(dataframe)
    n_records <- nrow(dt)

    min_year <- year(min(dt$start_date))
    max_year <- year(max(dt$end_date))
    cat(sprintf("Total records to process: %d\n", n_records))
    cat(sprintf("Date range in data: %d to %d\n", min_year, max_year))
    cat(sprintf("Batch size: %d records\n", batch_size))

    # --- File selection
    all_files <- list.files(path = input_folder_name, full.names = TRUE)
    pollutant_files <- all_files[grepl(paste0("^", pollutant), basename(all_files))]
    nc_files <- pollutant_files[grepl("\\.nc$", pollutant_files)]
    if (length(nc_files) == 0) stop(sprintf("No NetCDF files found for '%s'", pollutant))

    extract_years_from_filename <- function(filename) {
        year_range <- regmatches(filename, regexpr("\\d{4}-\\d{4}", filename))
        if (length(year_range) > 0) {
            years <- as.numeric(unlist(strsplit(year_range, "-")))
            return(list(start = years[1], end = years[2]))
        }
        single_years <- as.numeric(regmatches(filename, gregexpr("\\d{4}", filename))[[1]])
        valid_years <- single_years[single_years >= 1900 & single_years <= 2100]
        if (length(valid_years) > 0)
            return(list(start = min(valid_years), end = max(valid_years)))
        return(NULL)
    }
    # Find relevant files
    file_year_info <- lapply(nc_files, function(f) {
        yrs <- extract_years_from_filename(basename(f))
        overlaps <- !is.null(yrs) && !(yrs$end < min_year || yrs$start > max_year)
        list(file = f, overlaps = overlaps)
    })
    relevant_files <- nc_files[sapply(file_year_info, function(x) x$overlaps)]
    if (length(relevant_files) == 0) {
        relevant_files <- nc_files
        warning("No files matched year range; using all files.")
    }
    cat(sprintf("Loading %d relevant files...\n", length(relevant_files)))
    for (f in relevant_files) cat(sprintf("  - %s\n", basename(f)))

    cat("\nLoading raster data into memory...\n")
    raster_data <- rast(relevant_files)
    cat(sprintf("✓ Loaded %d layers\n", nlyr(raster_data)))

    all_dates <- as.Date(time(raster_data))
    if (any(is.na(all_dates))) stop("Found NA in raster dates")
    cat(sprintf("Date range in raster: %s to %s (%d days)\n",
        min(all_dates), max(all_dates), length(all_dates)))

    # Setup temp directory for intermediates
    if (is.null(temp_dir)) temp_dir <- file.path(output_folder_name, "temp")
    if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
    cat(sprintf("Temporary directory: %s\n", temp_dir))

    n_batches <- ceiling(n_records / batch_size)
    cat(sprintf("Processing %d batches of up to %d records each\n\n", n_batches, batch_size))
    all_temp_files <- character(n_batches)

    # Pre-compute numeric dates
    all_dates_num <- as.numeric(all_dates)

    for (batch_idx in seq_len(n_batches)) {
        batch_start_time <- Sys.time()
        batch_row_start <- (batch_idx - 1) * batch_size + 1
        batch_row_end <- min(batch_idx * batch_size, n_records)
        batch_rows <- batch_row_start:batch_row_end
        cat(sprintf("\n=== BATCH %d/%d (records %d-%d) ===\n",
            batch_idx, n_batches, batch_row_start, batch_row_end))

        batch_coords <- as.matrix(dt[batch_rows, .(long, lat)])
        cat("Extracting raster values...\n")
        extraction_start <- Sys.time()
        batch_extraction <- extract(raster_data, batch_coords)
        batch_matrix <- as.matrix(batch_extraction)
        extraction_time <- as.numeric(difftime(Sys.time(), extraction_start, units = "secs"))
        cat(sprintf("  Extraction completed in %.1f seconds\n", extraction_time))

        # Align matrix and dates
        if (ncol(batch_matrix) != length(all_dates)) {
            min_len <- min(ncol(batch_matrix), length(all_dates))
            batch_matrix <- batch_matrix[, 1:min_len, drop = FALSE]
            batch_dates_num <- all_dates_num[1:min_len]
        } else {
            batch_dates_num <- all_dates_num
        }

        batch_start_dates_num <- as.numeric(dt[batch_rows, start_date])
        batch_end_dates_num   <- as.numeric(dt[batch_rows, end_date])
        n_batch_records <- length(batch_rows)
        cat("Processing records in batch...\n")

        batch_date_masks <- sapply(seq_len(n_batch_records), function(i) {
            batch_dates_num >= batch_start_dates_num[i] & batch_dates_num <= batch_end_dates_num[i]
        }, simplify = "array")

        process_batch_record <- function(i) {
            mask <- batch_date_masks[, i]
            if (!any(mask)) return(numeric(0))
            vals <- batch_matrix[i, mask]
            vals <- vals[is.finite(vals)]
            if (length(vals) == 0) numeric(0) else round(vals, 2)
        }

        if (use_parallel && n_batch_records > 100) {
            batch_values <- parallel::mclapply(
                seq_len(n_batch_records),
                process_batch_record,
                mc.cores = n_cores,
                mc.preschedule = TRUE
            )
        } else {
            batch_values <- lapply(seq_len(n_batch_records), process_batch_record)
        }

        # Assemble batch result and save as RDS
        batch_result <- data.table(
            original_row = batch_rows,
            raster_value = batch_values,
            raster_n_values = sapply(batch_values, length),
            raster_has_data = sapply(batch_values, length) > 0
        )
        temp_file <- file.path(temp_dir, sprintf("batch_%04d.rds", batch_idx))
        saveRDS(batch_result, file = temp_file)
        all_temp_files[batch_idx] <- temp_file

        n_with_data <- sum(batch_result$raster_has_data)
        avg_values <- mean(batch_result$raster_n_values[batch_result$raster_has_data])
        batch_time <- as.numeric(difftime(Sys.time(), batch_start_time, units = "mins"))
        cat(sprintf("  Batch completed in %.1f minutes\n", batch_time))
        cat(sprintf("  Records with data: %d/%d (%.1f%%)\n", n_with_data, n_batch_records, 100 * n_with_data / n_batch_records))
        cat(sprintf("  Average values per record: %.1f\n", avg_values))
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        avg_per_batch <- elapsed / batch_idx
        eta <- avg_per_batch * (n_batches - batch_idx)
        cat(sprintf("  Total elapsed: %.1f min | ETA remaining: %.1f min\n", elapsed, eta))
        gc()
    }

    # Combine all batch results
    cat("\n=== COMBINING BATCH RESULTS ===\n")
    all_results <- lapply(all_temp_files, readRDS)
    combined_results <- rbindlist(all_results)
    setorder(combined_results, original_row) # restore order

    # Add to final dataframe
    dataframe$raster_value    <- combined_results$raster_value
    dataframe$raster_n_values <- combined_results$raster_n_values
    dataframe$raster_has_data <- combined_results$raster_has_data

    # Clean up temp
    cat("Cleaning up temporary files...\n")
    unlink(temp_dir, recursive = TRUE)

    # Save final result:
    if (!dir.exists(output_folder_name)) dir.create(output_folder_name, recursive = TRUE)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_filename <- sprintf("%s_extraction_%s.rds", pollutant, timestamp)
    output_path <- file.path(output_folder_name, output_filename)
    saveRDS(dataframe, file = output_path)

    end_time <- Sys.time()
    total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("LARGE-SCALE PROCESSING COMPLETE\n")
    cat(strrep("=", 70), "\n", sep = "")
    cat(sprintf("Output file: %s\n", output_path))
    cat(sprintf("Total records: %d\n", nrow(dataframe)))
    cat(sprintf("Records with data: %d (%.1f%%)\n",
        sum(dataframe$raster_has_data), 100 * mean(dataframe$raster_has_data)))
    cat(sprintf("Total processing time: %.1f hours (%.1f minutes)\n", total_time / 60, total_time))
    cat(sprintf("Average time per batch: %.1f minutes\n", total_time / n_batches))
    cat(sprintf("Estimated throughput: %.0f records/minute\n", nrow(dataframe) / total_time))
    cat(strrep("=", 70), "\n", sep = "")

    return(invisible(output_path))
}
