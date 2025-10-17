#' Function-level boxplots of pathway contributions
#'
#' Generate per-taxon boxplots of functional contributions (e.g., KO/EC/MetaCyc)
#' from stratified PICRUSt2 outputs, grouped by sample metadata categories.
#' This function merges ASV-level stratified contributions into a chosen
#' taxonomy level (e.g., Genus) and produces boxplots for functions that are
#' significant in differential abundance analysis results from
#' \code{pathway_sdaa()}.
#'
#' @param stratified_path Character scalar. Path to the stratified PICRUSt2 output
#'   (e.g., \code{pred_metagenome_contrib.tsv}), which must contain at least
#'   the columns \code{sample}, \code{taxon}, \code{function},
#'   \code{norm_taxon_function_contrib}.
#' @param pathway_sdaa_result Named list of results from \code{pathway_sdaa()}.
#'   Each element should be a data.frame with at least \code{feature} and
#'   \code{pathway_name} (or \code{feature} only). The list names must correspond
#'   to taxa labels (e.g., genera).
#' @param physeq A \code{phyloseq} object containing both taxonomy information
#'   (via \code{tax_table}) and sample metadata (via \code{sample_data}).
#' @param taxon_level Character scalar. Taxonomic rank at which to aggregate
#'   ASV-level contributions (e.g., \code{"Genus"}). Must be a column in
#'   \code{phyloseq::tax_table(physeq)}.
#' @param group Character scalar. Column name in \code{sample_data(physeq)}
#'   indicating the grouping variable (e.g., treatment, case/control).
#'
#'
#' @examples
#' \dontrun{
#' # Example usage
#' stratified_path <- "pred_metagenome_contrib.tsv"
#' plots <- function_box(
#'   stratified_path     = stratified_path,
#'   pathway_sdaa_result = sdaa_res,
#'   physeq              = physeq_obj,
#'   taxon_level         = "Genus",
#'   group               = "Group"
#' )
#' # Access the plot for GenusA
#' plots[["GenusA"]]
#' }
#'
#' @seealso \code{\link{pathway_sdaa}}, \pkg{phyloseq}, \pkg{ggplot2}
#'
#' @export
function_box <- function(
    stratified_path,
    pathway_sdaa_result,
    physeq,
    taxon_level = "Genus",
    group
){
  # --- checks -------------------------------------------------------------------
  stopifnot(is.character(stratified_path), length(stratified_path) == 1)
  if (!file.exists(stratified_path)) stop("File not found: ", stratified_path)
  stopifnot(inherits(physeq, "phyloseq"))
  if (!taxon_level %in% colnames(phyloseq::tax_table(physeq))) {
    stop("taxon_level '", taxon_level, "' not found in taxonomy table.")
  }
  if (!is.list(pathway_sdaa_result) || !length(pathway_sdaa_result)) {
    stop("`pathway_sdaa_result` must be a non-empty named list.")
  }
  if (is.null(names(pathway_sdaa_result)) || any(names(pathway_sdaa_result) == "")) {
    stop("`pathway_sdaa_result` must be a *named* list (names = genera).")
  }

  # --- IO: PICRUSt2 stratified table -------------------------------------------
  df <- data.table::fread(stratified_path, sep = "\t")
  req_cols <- c("sample","taxon","function","norm_taxon_function_contrib")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing columns in stratified file: ", paste(miss, collapse = ", "))

  # --- taxonomy & samples -------------------------------------------------------
  taxtab <- as.data.frame(phyloseq::tax_table(physeq))
  taxtab$ASV <- rownames(taxtab)
  taxtab <- taxtab[, c("ASV", taxon_level), drop = FALSE]

  sample_ids <- phyloseq::sample_names(physeq)
  if (!length(intersect(df$sample, sample_ids))) {
    stop("No overlapping samples between PICRUSt2 file and phyloseq object.")
  }

  meta <- as.data.frame(phyloseq::sample_data(physeq))
  if (!(group %in% colnames(meta))) {
    stop("Group column '", group, "' not found in sample_data.")
  }
  meta$sample <- rownames(meta)
  meta[[group]] <- as.factor(meta[[group]])

  # Helper to normalize function IDs (e.g., "ko:K00001" -> "K00001")
  .clean_fun <- function(x){
    x <- as.character(x)
    x <- sub("^ko:", "", x, ignore.case = TRUE)
    x <- sub("^KO:", "", x, ignore.case = TRUE)
    x
  }

  # --- ASV -> sample wide tables (rows = function IDs) --------------------------
  df_split <- split(df, df$taxon)

  per_asv <- lapply(df_split, function(dt){
    data.table::dcast(
      data.table::as.data.table(dt),
      `function` ~ sample,
      value.var = "norm_taxon_function_contrib",
      fill = 0
    )
  })

  # keep only columns present in physeq samples
  per_asv <- lapply(per_asv, function(dt){
    cols <- intersect(colnames(dt), c("function", sample_ids))
    dt[, ..cols]
  })

  # --- aggregate ASVs to chosen taxon_level ------------------------------------
  result_by_level <- list()
  for (asv in names(per_asv)) {
    lvl <- taxtab[match(asv, taxtab$ASV), taxon_level]
    if (length(lvl) == 0 || is.na(lvl) || lvl == "") next
    if (tolower(lvl) %in% c("unknown","uncultured")) next

    dt <- per_asv[[asv]]
    if (!"function" %in% colnames(dt) || ncol(dt) <= 1L) next

    # melt asv table
    dt_long <- data.table::melt(
      data.table::as.data.table(dt),
      id.vars = "function",
      variable.name = "sample",
      value.name = "value"
    )

    # append into level bucket
    if (is.null(result_by_level[[lvl]])) {
      result_by_level[[lvl]] <- dt_long
    } else {
      result_by_level[[lvl]] <- data.table::rbindlist(
        list(result_by_level[[lvl]], dt_long),
        use.names = TRUE, fill = TRUE
      )
    }
  }

  # sum duplicates (function x sample) within each level and recast to wide
  for (lvl in names(result_by_level)) {
    long <- result_by_level[[lvl]]
    if (!nrow(long)) { result_by_level[[lvl]] <- NULL; next }
    agg <- long[, .(value = sum(value, na.rm = TRUE)), by = .(`function`, sample)]
    wide <- data.table::dcast(agg, `function` ~ sample, value.var = "value", fill = 0)
    result_by_level[[lvl]] <- wide
  }

  # --- build plots by genus present in pathway_sdaa_result ----------------------
  genera <- names(pathway_sdaa_result)
  plot_list <- vector("list", length(genera)); names(plot_list) <- genera

  for (genus in genera) {
    genus_tbl <- result_by_level[[genus]]
    if (is.null(genus_tbl) || !"function" %in% colnames(genus_tbl) || ncol(genus_tbl) <= 1L) {
      next
    }

    # long format (function × sample)
    g_long <- data.table::melt(
      data.table::as.data.table(genus_tbl),
      id.vars = "function", variable.name = "sample", value.name = "abundance"
    )
    if (!nrow(g_long)) next

    # normalize function IDs to match DAA result "feature"
    g_long[, function_clean := .clean_fun(`function`)]

    # DAA table for this genus
    daa_df <- pathway_sdaa_result[[genus]]
    if (is.null(daa_df) || !nrow(daa_df)) next

    # prefer pathway_name, fall back to feature
    if (!"feature" %in% names(daa_df)) next
    if (!"pathway_name" %in% names(daa_df)) daa_df$pathway_name <- daa_df$feature

    daa_df$feature_clean <- .clean_fun(daa_df$feature)
    daa_keep <- daa_df[, c("feature_clean","pathway_name"), drop = FALSE]

    # join to keep only functions present in DAA result and carry names
    g_long <- merge(
      g_long,
      daa_keep,
      by.x = "function_clean",
      by.y = "feature_clean",
      all.x = FALSE, all.y = FALSE
    )
    if (!nrow(g_long)) next

    # attach metadata and filter to samples with group
    g_long <- merge(g_long, meta, by = "sample", all.x = TRUE)
    g_long <- g_long[!is.na(g_long[[group]]), , drop = FALSE]
    if (!nrow(g_long)) next

    # optional: log-scale helper if you need it later (kept numeric safe)
    # eps <- 1e-6
    # g_long$abundance_log <- log10(g_long$abundance + eps)

    # plot
    p <- ggplot2::ggplot(
      g_long,
      ggplot2::aes(x = stats::reorder(pathway_name, abundance, FUN = median),  # order by median
                   y = abundance, fill = .data[[group]])
    ) +
      ggplot2::geom_boxplot(outlier.shape = 16, outlier.size = 0.8) +
      ggplot2::labs(
        title = genus,
        x = "Function",
        y = "Contribution",
        fill = group
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title.position = "plot",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    plot_list[[genus]] <- p
  }

  # drop NULLs (genera with no plot)
  plot_list[!vapply(plot_list, is.null, logical(1))]
}



