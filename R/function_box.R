#' Plot Functional Abundance Boxplots for Driving Taxa
#'
#' This function generates boxplots of stratified functional pathway abundances
#' (from PICRUSt2 output) at a specified taxonomic level for each driving taxon,
#' grouped by a metadata variable.
#'
#' @param stratified_path Path to the stratified output file from PICRUSt2 (typically `*_strat.tsv`).
#' @param pathway_sdaa_result A named list of data.frames containing differential results
#'   (with `feature` and `pathway_name`) for each taxon.
#' @param physeq A `phyloseq` object containing the microbial abundance table, taxonomy, and sample metadata.
#' @param taxon_level A string indicating the taxonomic level to aggregate functions by
#'   (e.g., `"Genus"`, `"Family"`). Default is `"Genus"`.
#' @param group A string specifying the sample metadata column used for grouping in boxplots.
#'
#' @return A named list of `ggplot2` boxplot objects, one for each taxon with matching functional results.
#'
#' @details This function:
#' \itemize{
#'   \item Parses PICRUSt2 stratified output and maps ASVs to their taxonomic levels.
#'   \item Aggregates function abundances at the taxon level.
#'   \item Filters and matches functions found in the user-provided differential result per taxon.
#'   \item Generates and returns a list of ggplot boxplots of function abundance by group.
#' }
#'
#' @import ggplot2 phyloseq data.table reshape2
#' @importFrom utils read.table
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' \dontrun{
#' plot_list <- function_box(
#'   stratified_path = "ko_metagenome_out/pred_metagenome_unstrat.tsv",
#'   pathway_sdaa_result = diff_results,
#'   physeq = my_phyloseq,
#'   taxon_level = "Genus",
#'   group = "disease_status"
#' )
#' plot_list$Bacteroides  # View boxplot for Bacteroides genus
#' }
function_box <- function(stratified_path, pathway_sdaa_result, physeq, taxon_level = "Genus", group) {
  # library(ggplot2)
  # library(phyloseq)
  # library(data.table)
  # library(reshape2)

  # Load stratified Picrust2 output
  df <- fread(stratified_path, sep = "\t")

  # Merge taxonomy to desired level
  tax_table_df <- as.data.frame(tax_table(physeq))
  tax_table_df$ASV <- rownames(tax_table_df)
  tax_table_df <- tax_table_df[, c("ASV", taxon_level)]

  df_split <- split(df, df$taxon)
  result_list <- lapply(df_split, function(dt) {
    data.table::dcast(dt, `function` ~ sample, value.var = "taxon_rel_function_abun", fill = 0)
  })

  sample_names <- sample_names(physeq)
  result_list_subset <- lapply(result_list, function(dt) {
    cols_to_keep <- c("function", intersect(names(dt), sample_names))
    dt[, ..cols_to_keep]
  })

  result_list_merged <- list()
  for (i in seq_along(result_list_subset)) {
    asv <- names(result_list_subset)[i]
    level <- tax_table_df[tax_table_df$ASV == asv, taxon_level]
    if (length(level) == 0 || is.na(level) || level == "") {
      level <- "unknown"
    }
    dt <- result_list_subset[[i]]

    if (!"function" %in% colnames(dt)) next

    dt_long <- melt(dt, id.vars = "function")
    colnames(dt_long) <- c("function", "sample", "value")

    if (level %in% names(result_list_merged)) {
      existing <- result_list_merged[[level]]
      existing_long <- melt(existing, id.vars = "function")
      colnames(existing_long) <- c("function", "sample", "value")

      combined <- rbind(existing_long, dt_long)
      combined <- aggregate(value ~ `function` + sample, data = combined, sum)
      result_list_merged[[level]] <- dcast(combined, `function` ~ sample, value.var = "value", fill = 0)
    } else {
      result_list_merged[[level]] <- dt
    }
  }

  meta <- data.frame(sample_data(physeq))
  meta$sample <- rownames(meta)

  plot_list <- list()
  genera <- names(pathway_sdaa_result)

  for (genus in genera) {
    genus_df <- result_list_merged[[genus]]
    genus_diff <- pathway_sdaa_result[[genus]]

    # Prepare long-format data for all functions in the genus
    genus_df_long <- melt(genus_df, id.vars = "function", variable.name = "sample", value.name = "abundance")
    genus_df_long$function_clean <- gsub("^ko:", "", genus_df_long$`function`)
    genus_df_long <- merge(genus_df_long, meta, by = "sample")

    # Only plot functions that exist in the differential result
    # Extract numeric parts
    genus_df_long$function_num <- gsub("^[A-Za-z]+", "", genus_df_long$function_clean)
    genus_diff$feature_num <- gsub("^ko", "", genus_diff$feature)

    # Filter based on numeric part match
    # Merge pathway_name info into genus_df_long
    genus_df_long <- merge(genus_df_long, genus_diff[, c("feature_num", "pathway_name")],
                           by.x = "function_num", by.y = "feature_num", all.x = TRUE)
    # Remove NA or empty pathway names
    genus_df_long <- genus_df_long[!is.na(genus_df_long$pathway_name) & genus_df_long$pathway_name != "", ]

    if (nrow(genus_df_long) == 0) next


    p <- ggplot(genus_df_long, aes(x = pathway_name, y = abundance, fill = get(group))) +
      geom_boxplot() +
      labs(title = genus, x = "Function", y = "Abundance",fill=group) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list[[genus]] <- p
  }

  return(plot_list)
}


