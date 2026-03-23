naive_mean <- function(df, domain_var, y) {
  
  # Remove rows where y is NA
  df <- df[!is.na(df[[y]]), , drop = FALSE]
  y_vals <- df[[y]]
  
  # If no domain variable is provided, return overall mean and SE
  if (is.null(domain_var)) {
    n <- length(y_vals)
    return(list(
      mean   = if (n > 0) mean(y_vals) else NA_real_,
      var    = if (n > 1) var(y_vals) / n else NA_real_
    ))
  }
  
  # Check if domain_var exists in the data frame
  if (!domain_var %in% names(df)) {
    stop("domain_var '", domain_var, "' not found in df.")
  }
  
  # Check if domain_var contains NA values
  if (any(is.na(df[[domain_var]]))) {
    stop("NA found in domain_var column: ", domain_var)
  }
  
  # Calculate mean and SE for rows where domain_var == 1
  idx <- df[[domain_var]] == 1
  domain_data <- y_vals[idx]
  n <- length(domain_data)
  
  return(list(mean   = if (n > 0) mean(domain_data) else NA_real_, 
              var    = if (n > 1) var(domain_data) / n else NA_real_))
  
}

