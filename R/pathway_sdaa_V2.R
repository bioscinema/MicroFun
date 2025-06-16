#' Perform Stratified Differential Abundance Analysis of Functional Pathways
#'
#' This function performs stratified differential abundance analysis (DAA) using
#' PICRUSt2 functional output and phyloseq metadata, either by grouping ASV contributions
#' to functions (direction = "taxa") or by grouping taxonomic contributions to a specific function (direction = "function").
#'
#' @param stratified_path Path to the stratified output file from PICRUSt2 (typically `*_strat.tsv`).
#' @param physeq A `phyloseq` object containing sample data and taxonomic information.
#' @param taxon_level Character. Taxonomic level to aggregate contributions (e.g., `"Genus"`, `"Family"`). Default is `"Genus"`.
#' @param group Character. Column name in `sample_data(physeq)` that defines the grouping variable.
#' @param daa_method Character. Differential abundance method to use, e.g., `"Maaslin2"`, `"ALDEx2"`, etc. Default is `"Maaslin2"`.
#' @param p.adjust.method Character. Method for p-value adjustment (e.g., `"fdr"`, `"none"`). Default is `"none"`.
#' @param threshold Numeric. Significance threshold for adjusted p-values. Default is `0.05`.
#' @param direction Character. Whether to group by `"taxa"` or `"function"` for stratification. Default is `"taxa"`.
#'
#' @return A named list of annotated differential abundance results per taxon or function.
#'
#' @details
#' When `direction = "taxa"`, the function:
#' \itemize{
#'   \item Aggregates KO/function contributions from each taxon.
#'   \item Maps ASVs to taxonomic levels.
#'   \item Performs DAA using the selected method.
#'   \item Annotates significant features with KEGG pathway names.
#' }
#' When `direction = "function"`, the function:
#' \itemize{
#'   \item Aggregates taxon contributions per function.
#'   \item Compares taxa for each function using DAA.
#'   \item Annotates results using pathway names if available.
#' }
#'
#' @import data.table
#' @importFrom reshape2 melt dcast
#' @importFrom phyloseq sample_names sample_data tax_table
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pathway_sdaa(
#'   stratified_path = "pred_metagenome_strat.tsv",
#'   physeq = my_phyloseq_object,
#'   taxon_level = "Genus",
#'   group = "disease_status",
#'   daa_method = "Maaslin2",
#'   p.adjust.method = "fdr",
#'   threshold = 0.05,
#'   direction = "taxa"
#' )
#' }
pathway_sdaa_v2 <- function(
    stratified_path,
    physeq,
    taxon_level = "Genus",
    group,
    daa_method = "Maaslin2",
    p.adjust.method = "none",
    threshold=0.05,
    direction = "taxa"
) {
  # Load stratified Picrust2 output
  df <- data.table::fread(stratified_path, sep = "\t")

  sample_names <- phyloseq::sample_names(physeq)
  mysam <- data.frame(physeq@sam_data)

  if (direction == "taxa") {

    # Split by taxon
    taxon_list <- split(df, df$taxon)

    # Convert each taxon into a wide-format table (KO x samples)
    result_list <- lapply(taxon_list, function(dt) {
      data.table::dcast(dt, `function` ~ sample, value.var = "taxon_function_abun", fill = 0)
    })
    # names(result_list) <- names(taxon_list)

    # Filter to samples that exist in physeq

    result_list_subset <- lapply(result_list, function(dt) {
      # print(names(dt))
      cols_to_keep <- c("function", intersect(names(dt), sample_names))
      dt_subset <- dt[, ..cols_to_keep]  # data.table column subset
      return(dt_subset)
    })

    result_list_subset <- result_list_subset[
      sapply(result_list_subset, ncol) > 1
    ]

    # Merge taxonomy to desired level if specified
    if (!is.null(taxon_level)) {
      tax_table_df <- as.data.frame(phyloseq::tax_table(physeq))
      tax_table_df$ASV <- rownames(tax_table_df)
      tax_table_df <- tax_table_df[, c("ASV", taxon_level)]

      result_list_merged <- list()

      for (i in seq_along(result_list_subset)) {
        asv <- names(result_list_subset)[i]
        level <- tax_table_df[tax_table_df$ASV == asv, taxon_level]
        if (length(level) == 0 || is.na(level) || level == "" ) {
          next
        }
        if (tolower(level) == "unknown" || tolower(level) == "uncultured") next
        # level <- ifelse(is.na(level) | length(level) == 0, asv, level)
        # print(level)
        dt <- result_list_subset[[i]]

        if (!"function" %in% colnames(dt)) next

        dt_long <- as.data.frame(reshape2::melt(dt, id.vars = "function"))
        colnames(dt_long) <- c("function", "sample", "value")

        if (level %in% names(result_list_merged)) {
          existing <- result_list_merged[[level]]
          if (!"function" %in% colnames(existing)) next

          existing_long <- as.data.frame(reshape2::melt(existing, id.vars = "function"))
          colnames(existing_long) <- c("function", "sample", "value")

          combined_long <- rbind(existing_long, dt_long)
          combined_long <- aggregate(value ~ `function` + sample, data = combined_long, sum)

          merged_dt <- reshape2::dcast(combined_long, `function` ~ sample, value.var = "value", fill = 0)
          result_list_merged[[level]] <- merged_dt
        } else {
          result_list_merged[[level]] <- dt
        }
      }

      result_list_merged <- lapply(result_list_merged, function(tbl) {
        tbl[is.na(tbl)] <- 0
        return(tbl)
      })

      # result_list_subset <- result_list_merged
      # Filter to retain only taxa with at least two groups represented

      result_list_subset <- list()

      for (genus in names(result_list_merged)) {
        tbl        <- result_list_merged[[genus]]
        # get your sample columns, excluding the first "function" column
        sample_ids <- setdiff(colnames(tbl), "function")
        # drop any samples with missing metadata
        sample_ids <- intersect(sample_ids, rownames(mysam))

        # look up their group labels
        sample_groups <- mysam[sample_ids, group]
        # drop NAs
        valid        <- !is.na(sample_groups)
        sample_ids   <- sample_ids[valid]
        sample_groups<- sample_groups[valid]

        # check #1: need at least two distinct groups
        if (length(unique(sample_groups)) < 2) {
          message("Skipping ", genus, "-only one group present")
          next
        }

        # check #2: need more than 4 total samples
        if (length(sample_ids) <= 4) {
          message("Skipping ", genus,
                  " -only ", length(sample_ids), " samples (need >4)")
          next
        }

        # if we reach here, both checks passed
        result_list_subset[[genus]] <- tbl
      }

      # result_list_subset <- lapply(result_list_subset, function(tbl) {
      #   # keep the "function" column aside
      #   fun <- tbl[, 1, drop = FALSE]
      #   mat <- as.matrix(tbl[, -1, drop = FALSE])
      #
      #   # compute column sums
      #   cs  <- colSums(mat, na.rm = TRUE)
      #   cs[cs == 0] <- NA  # avoid divide-by-zero
      #
      #   # divide each column by its sum
      #   mat_rel <- sweep(mat, 2, cs, `/`)
      #   mat_rel[is.na(mat_rel)] <- 0  # turn NA back to 0
      #
      #   # reassemble
      #   out <- data.frame(fun, mat_rel, check.names = FALSE, row.names = NULL)
      #   names(out)[1] <- "function"
      #   return(out)
      # })

    }
  } else if (direction == "function") {

    # Group by function
    function_list <- split(df, df$'function')
    result_list_subset <- list()

    # Extract taxonomy table
    tax_table_df <- as.data.frame(phyloseq::tax_table(physeq))
    tax_table_df$ASV <- rownames(tax_table_df)

    # Check taxon level exists
    if (!(taxon_level %in% colnames(tax_table_df))) {
      stop("The specified 'taxon_level' is not found in the taxonomy table.")
    }

    for (func_name in names(function_list)) {
      dt <- function_list[[func_name]]

      # Remove unclassified entries
      dt <- dt[dt$taxon != "unclassified", ]

      # Map ASV/taxon to target level
      dt$level <- tax_table_df[match(dt$taxon, tax_table_df$ASV), taxon_level]

      # Filter low-quality entries
      dt <- dt[!is.na(dt$level) & dt$level != "" &
                 !(tolower(dt$level) %in% c("unknown", "uncultured")), ]

      # Aggregate by taxon_level and sample
      dt_agg <- dt[, .(norm_taxon_function_contrib = sum(norm_taxon_function_contrib, na.rm = TRUE)),
                   by = .(level, sample)]

      # Convert to wide format: rows = genus, cols = samples
      dt_wide <- data.table::dcast(dt_agg, level ~ sample,
                                   value.var = "norm_taxon_function_contrib", fill = 0)

      # Keep only sample columns in both dt_wide and metadata
      sample_ids <- intersect(colnames(dt_wide), sample_names)
      dt_wide <- dt_wide[, c("level", sample_ids), with = FALSE]

      # Filter groups
      present_samples <- intersect(sample_ids, rownames(mysam))
      sample_groups <- mysam[present_samples, group]
      valid <- !is.na(sample_groups)
      sample_groups <- sample_groups[valid]
      sample_ids <- present_samples[valid]

      if (length(unique(sample_groups)) < 2) next
      if (length(sample_ids) <= 4) next

      # Subset again to filtered samples
      dt_wide <- dt_wide[, c("level", sample_ids), with = FALSE]

      # Normalize by sample (TSS)
      # taxon <- dt_wide[, 1, drop = FALSE]
      # mat <- as.matrix(dt_wide[, -1, with = FALSE])
      # cs <- colSums(mat, na.rm = TRUE)
      # cs[cs == 0] <- NA
      # mat_rel <- sweep(mat, 2, cs, `/`)
      # mat_rel[is.na(mat_rel)] <- 0

      result_tbl <- data.frame(dt_wide)
      names(result_tbl)[1] <- taxon_level

      result_list_subset[[func_name]] <- result_tbl
    }
  }




  # Initialize result list
  all_daa_results <- list()

  # Iterate through taxon-level KO tables with a progress bar
  total <- length(result_list_subset)
  for (i in seq_along(result_list_subset)) {
    name <- names(result_list_subset)[i]
    cat(sprintf("\n[%d/%d] Processing %s: %s\n", i, total,
                ifelse(direction == "taxa", "taxon", "function"), name))

    cat("Running differential analysis: [", sep = "")
    for (j in 1:30) {
      cat("=")
      flush.console()
      Sys.sleep(0.01)  # simulate progress bar animation
    }
    cat("]\n")

    df1 <- data.frame(result_list_subset[[name]], check.names = FALSE)

    if (direction == "taxa") {
      # Remove "ko:" prefix from function names
      df1$`function` <- gsub("^ko:", "", df1$`function`)

      # Convert from KO-level to KEGG pathway abundance
      ko_abun <- ggpicrust2::ko2kegg_abundance(data = df1)

      # Prepare metadata
      mysam_subset <- mysam[rownames(mysam) %in% colnames(ko_abun), ]
      mysam_subset$sample <- rownames(mysam_subset)
      mysam_subset[[group]] <- as.factor(mysam_subset[[group]])

      # DAA
      result <- pathway_daa(
        abundance = ko_abun,
        metadata = mysam_subset,
        group = group,
        daa_method = daa_method,
        p.adjust = p.adjust.method
      )

      # Annotate and store result
      annotated <- tryCatch(
        {
          pathway_annotation(
            pathway        = "KO",
            daa_results_df = result,
            ko_to_kegg     = TRUE,
            threshold      = threshold
          )
        },
        error = function(e) {
          message("  -> Skipping ", name, ": ", e$message)
          NULL
        }
      )

    } else if (direction == "function") {
      # Remove "ko:" from row name (taxon column contains the taxon, rows are KO-specific already)
      # df1$taxon <- gsub("^ko:", "", df1$taxon)
      rownames(df1) <- df1[,taxon_level]
      df1 <- df1[, -which(colnames(df1) == taxon_level)]

      # Prepare metadata
      mysam_subset <- mysam[rownames(mysam) %in% colnames(df1), ]
      mysam_subset$sample <- rownames(mysam_subset)
      mysam_subset[[group]] <- as.factor(mysam_subset[[group]])

      # DAA directly on function-level abundance
      result <- pathway_daa(
        abundance = df1,
        metadata = mysam_subset,
        group = group,
        daa_method = daa_method,
        p.adjust = p.adjust.method
      )

      annotated <- tryCatch(
        {
          res <- result_annotation(result, threshold = threshold)
          if (is.null(res) || nrow(res) == 0) {
            message("  -> Skipping ", name, ": no significant features found")
            res <- NULL
          }
          res
        },
        error = function(e) {
          message("  -> Skipping ", name, ": ", e$message)
          NULL
        }
      )

    }

    # 5) Store if non-NULL
    if (!is.null(annotated)) {
      all_daa_results[[name]] <- annotated
    }
  }
  if (length(all_daa_results) == 0) {
    stop(sprintf("No statistically significant biomarkers found for any %s.\nTry a less stringent threshold or review your input data.",
                 ifelse(direction == "taxa", "taxon", "function")))
  }
  return(all_daa_results)
}


