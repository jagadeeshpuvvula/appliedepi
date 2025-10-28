getRasterValueInRange_Fast <- function(dataframe, input_folder_name, pollutant = "PM") {
    dt <- as.data.table(dataframe)
    
    # Get all NetCDF files that START with the pollutant name
    all_files <- list.files(path = input_folder_name, full.names = TRUE)
    pollutant_files <- all_files[grepl(paste0("^", pollutant), basename(all_files))]
    nc_files <- pollutant_files[grepl("\\.nc$", pollutant_files)]
    
    if (length(nc_files) == 0) {
        stop(paste("No NetCDF files found starting with:", pollutant))
    }
    
    print(paste("Found", length(nc_files), "NetCDF files"))
    
    # Load all NetCDF files
    raster_data <- rast(nc_files)
    all_dates <- time(raster_data)
    
    print(paste("Loaded", nlyr(raster_data), "layers"))
    
    # Extract values for ALL locations at once (FASTEST approach)
    coordinates <- as.matrix(dt[, .(long, lat)])
    
    print("Extracting values for all locations...")
    all_extractions <- extract(raster_data, coordinates)
    
    print("Processing daily values...")
    
    # Convert to regular dataframe to avoid data.table scoping issues in parallel
    df <- as.data.frame(dataframe)
    
    # Setup parallel processing
    num_cores <- min(detectCores() - 7, 25)
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    clusterExport(cl, varlist = c("all_extractions", "all_dates", "df"), 
                  envir = environment())
    
    # Process each row
    raster_values <- foreach(i = 1:nrow(df), 
                             .packages = c("lubridate")) %dopar% {
                                 start_date <- as.Date(df$start_date[i])
                                 end_date <- as.Date(df$end_date[i])
                                 
                                 # Get extracted values for this location
                                 location_values <- as.numeric(all_extractions[i, -1])  # Remove ID column
                                 
                                 # Find indices for dates in the range
                                 date_idx <- which(all_dates >= start_date & all_dates <= end_date)
                                 
                                 if (length(date_idx) == 0) {
                                     return(NULL)
                                 }
                                 
                                 # Get daily values for the date range
                                 daily_values <- location_values[date_idx]
                                 
                                 # Clean the values (remove NA, NaN, Inf)
                                 daily_values <- daily_values[!is.na(daily_values) & 
                                                              !is.nan(daily_values) & 
                                                              is.finite(daily_values)]
                                 
                                 if (length(daily_values) == 0) {
                                     return(NULL)
                                 }
                                 
                                 # Round to 2 decimal places
                                 daily_values <- round(daily_values, 2)
                                 
                                 return(daily_values)
                             }
    
    stopCluster(cl)
    
    dataframe$raster_value <- raster_values
    
    return(dataframe)
}
