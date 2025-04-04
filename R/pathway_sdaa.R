#' @title Stratified Pathway Differential Abundance Analysis
#' @description Performs stratified differential abundance analysis (sDAA) using stratified KO predictions from PICRUSt2.
#' For each taxon, this function extracts its KO function profile and applies a differential abundance method (e.g., MaAsLin2).
#' Optionally, taxa can be aggregated to a specified taxonomic level before analysis.
#'
#' @param picrust2_path Path to PICRUSt2 stratified KO output file ("pred_metagenome_contrib.tsv").
#' @param physeq A phyloseq object containing sample metadata and taxonomy.
#' @param taxon_level Optional. Taxonomic level to aggregate taxa (e.g., "Genus"). If NULL, no aggregation is done.
#' @param group Character. The column name in metadata to use for group comparison (e.g., "prediabetes").
#' @param daa_method Character. Differential abundance method. Default is "Maaslin2".
#' @param p.adjust.method Character. Multiple testing correction method. Default is "BH".
#' @param output_dir Character. Directory to store MaAsLin2 output. Default is "maaslin2_output".
#'
#' @return A named list of data frames containing KEGG-annotated significant features for each taxon.
#' @export
pathway_sdaa <- function(
    stratified_path,
    physeq,
    taxon_level = NULL,
    group,
    daa_method = "Maaslin2",
    p.adjust.method = "BH"
) {
  # Load stratified Picrust2 output
  df <- data.table::fread(stratified_path, sep = "\t")

  # Split by taxon
  taxon_list <- split(df, df$taxon)

  # Convert each taxon into a wide-format table (KO x samples)
  result_list <- lapply(taxon_list, function(dt) {
    data.table::dcast(dt, `function` ~ sample, value.var = "taxon_rel_function_abun", fill = 0)
  })
  names(result_list) <- names(taxon_list)

  # Filter to samples that exist in physeq
  sample_names <- phyloseq::sample_names(physeq)
  result_list_subset <- lapply(result_list, function(dt) {
    cols_to_keep <- c("function", intersect(names(dt), sample_names))
    dt_subset <- dt[, ..cols_to_keep]  # data.table column subset
    return(dt_subset)
  })

  # Merge taxonomy to desired level if specified
  if (!is.null(taxon_level)) {
    tax_table_df <- as.data.frame(phyloseq::tax_table(physeq))
    tax_table_df$ASV <- rownames(tax_table_df)
    tax_table_df <- tax_table_df[, c("ASV", taxon_level)]

    merged_names <- sapply(names(result_list_subset), function(x) {
      matched_level <- tax_table_df[tax_table_df$ASV == x, taxon_level]
      if (length(matched_level) == 0 || is.na(matched_level)) {
        return(x)
      } else {
        return(matched_level)
      }
    })
    names(result_list_subset) <- merged_names
  }

  # Initialize result list
  all_daa_results <- list()

  # Iterate through taxon-level KO tables with a progress bar
  total <- length(result_list_subset)
  for (i in seq_along(result_list_subset)) {
    name <- names(result_list_subset)[i]
    cat(sprintf("\n[%d/%d] Processing taxon: %s\n", i, total, name))
    cat("Running differential analysis: [", sep = "")
    for (j in 1:30) {
      cat("=")
      flush.console()
      Sys.sleep(0.01)  # simulate progress bar animation
    }
    cat("]\n")

    df1 <- data.frame(result_list_subset[[name]], check.names = FALSE)
    df1$`function` <- gsub("^ko:", "", df1$`function`)
    ko_abun <- ggpicrust2::ko2kegg_abundance(data = df1)

    # Match sample metadata
    mysam <- data.frame(physeq@sam_data)
    mysam_subset <- mysam[rownames(mysam) %in% colnames(df1), ]

    # Run differential analysis using pathway_daa
    result <- ggpicrust2::pathway_daa(
      abundance = ko_abun,
      metadata = mysam_subset,
      group = group,
      daa_method = daa_method,
      p.adjust = p.adjust.method
    )

    # Annotate and store result
    annotated <- ggpicrust2::pathway_annotation("KO", result, ko_to_kegg = TRUE)
    all_daa_results[[name]] <- annotated
  }

  return(all_daa_results)
}




