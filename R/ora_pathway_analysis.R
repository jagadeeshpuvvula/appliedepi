#usage
#results <- ora_pathway_enrichment(res_df,
#                                  metabolite_fdr_col = "fdr",
#                                  metabolite_fdr_threshold = 0.2,
#                                  pathway_col = "sub_pathway",
#                                  min_pathway_size = 2,
#                                  min_sig_in_pathway = 1,
#                                  stratify_by_sex = TRUE)


ora_pathway_enrichment <- function(res_df, 
                                            metabolite_fdr_col = "fdr",
                                            metabolite_fdr_threshold = 0.2,
                                            pathway_col = "sub_pathway",
                                            min_pathway_size = 2,
                                            min_sig_in_pathway = 1,
                                            stratify_by_sex = TRUE) {
  
  # Validate inputs
  if(metabolite_fdr_threshold < 0 || metabolite_fdr_threshold > 1){
    stop("metabolite_fdr_threshold must be between 0 and 1")
  }
  
  # Check required columns
  required_cols <- c("variable", "chem_id", "chemical_name", "estimate", metabolite_fdr_col, pathway_col)
  if(stratify_by_sex) {
    required_cols <- c(required_cols, "sex")
  }
  
  missing_cols <- setdiff(required_cols, colnames(res_df))
  if(length(missing_cols) > 0){
    stop("Input data frame missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get all variable-(sex) combinations
  if(stratify_by_sex) {
    combinations <- res_df %>%
      distinct(variable, sex) %>%
      arrange(variable, sex)
  } else {
    combinations <- res_df %>%
      distinct(variable) %>%
      mutate(sex = "all") %>%
      arrange(variable)
  }
  
  single_enrichment <- function(var_name, sex_name) {
    # Filter data based on stratification
    if(stratify_by_sex) {
      var_data <- res_df %>% filter(variable == var_name, sex == sex_name)
    } else {
      var_data <- res_df %>% filter(variable == var_name)
    }
    
    # Background metabolites for combination
    background <- var_data %>%
      distinct(chem_id)
    M <- nrow(background)
    
    # Significant metabolites in background (by dynamic FDR column)
    sig_background <- var_data %>%
      filter(!!sym(metabolite_fdr_col) < metabolite_fdr_threshold) %>%
      distinct(chem_id)
    n <- nrow(sig_background)
    
    # Skip if no significant metabolites
    if (n == 0) {
      return(tibble())
    }
    
    # Pathway metabolites for combination
    pathway_data <- var_data %>%
      distinct(chem_id, chemical_name, estimate, !!sym(pathway_col)) %>%
      filter(!is.na(!!sym(pathway_col)))
    
    all_pathways <- pathway_data %>%
      pull(!!sym(pathway_col)) %>%
      unique()
    
    # Calculate enrichment for each pathway
    result <- map_df(all_pathways, function(pathway) {
      pathway_mets <- pathway_data %>%
        filter(!!sym(pathway_col) == pathway) %>%
        pull(chem_id)
      
      N <- length(pathway_mets)
      
      sig_pathway_mets <- intersect(sig_background$chem_id, pathway_mets)
      k <- length(sig_pathway_mets)
      
      expected <- n * (N / M)
      enrichment_ratio <- ifelse(expected > 0, k / expected, NA)
      
      # Calculate p-value using hypergeometric test
      p_value <- if(k > 0 && N > 0 && N <= M && n <= M){
        phyper(q = k - 1, m = N, n = M - N, k = n, lower.tail = FALSE)
      } else {
        1
      }
      
      # Get information about significant metabolites
      sig_met_info <- var_data %>%
        filter(chem_id %in% sig_pathway_mets) %>%
        arrange(!!sym(metabolite_fdr_col)) %>%
        mutate(met_fdr = paste0(chemical_name, " (est=", round(estimate, 3), 
                                ", FDR=", format(!!sym(metabolite_fdr_col), digits = 2), ")")) %>%
        pull(met_fdr) %>%
        paste(collapse = "; ")
      
      tibble(
        variable = var_name,
        sex = sex_name,
        pathway = pathway,
        total_in_pathway = N,
        significant_in_pathway = k,
        total_tested = M,
        total_significant = n,
        expected = expected,
        enrichment_ratio = enrichment_ratio,
        p_value = p_value,
        significant_metabolites = sig_met_info
      )
    }) %>%
      filter(total_in_pathway >= min_pathway_size, significant_in_pathway >= min_sig_in_pathway) %>%
      mutate(
        fdr = p.adjust(p_value, method = "fdr")
      ) %>%
      arrange(p_value)
    
    return(result)
  }
  
  # Perform analysis for all combinations
  all_results <- map2_df(combinations$variable, combinations$sex, ~single_enrichment(.x, .y))
  
  return(all_results)
}
