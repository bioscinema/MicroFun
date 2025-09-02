#' Annotate Significant Features from Differential Abundance Results from function direction
#'
#' This function filters and returns statistically significant features from a differential 
#' abundance analysis (DAA) result based on an adjusted p-value threshold.
#'
#' @param daa_result A data frame containing the results of a differential abundance analysis. 
#'   Must include at least the columns \code{feature} and \code{p_adjust}.
#' @param threshold A numeric value specifying the adjusted p-value cutoff for significance. 
#'   Defaults to \code{0.05}.
#'
#' @return A data frame of significant features with adjusted p-values below the given threshold, 
#'   sorted by ascending p-value.
#'
#'
#' @keywords internal
result_annotation <- function(daa_result, threshold = 0.05) {
  # Check required columns
  required_cols <- c("feature", "p_adjust")
  if (!all(required_cols %in% colnames(daa_result))) {
    stop("Input must contain at least 'feature' and 'p_adjust' columns.")
  }
  
  # Subset to significant features
  significant <- daa_result[!is.na(daa_result$p_adjust) & daa_result$p_adjust < threshold, ]
  
  # Optional: sort by adjusted p-value
  significant <- significant[order(significant$p_adjust), ]
  
  return(significant)
}

