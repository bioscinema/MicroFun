#' Plot function-to-taxon contributions as grouped boxplots
#'
#' This function takes a stratified PICRUSt2 output table, a phyloseq object,
#' and a list of pathway-level results, and produces ggplot2 boxplots of
#' function contributions summarized by taxonomic level and experimental group.
#'
#' @param stratified_path Character string. Path to a stratified PICRUSt2 output file
#'   (typically `*_stratified.tsv`), containing columns for function, taxon, sample, and contribution values.
#' @param pathway_sdaa_result A named list of results from an SDAA (statistical differential abundance analysis)
#'   step, keyed by function IDs. Only names are used here.
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object containing sample metadata
#'   and taxonomy table. Sample IDs must match those in \code{stratified_path}.
#' @param group Character scalar. Column name in sample metadata to use for grouping (e.g. treatment).
#' @param taxon_level Character scalar. It would be better to align your taxonomy level
#' to your pathway_sdaa_result.
#' @param pathway Character scalar. Functional annotation scheme used: one of
#'   \code{"KO"}, \code{"EC"}, or \code{"MetaCyc"}. Defaults to "KO".
#' @param value_col Character scalar. Name of the numeric contribution column
#'   in the stratified table. Defaults to \code{"norm_taxon_function_contrib"}.
#' @param wrap_width Integer. Width (in characters) for line-wrapping taxon labels
#'   on the x-axis. Defaults to 36.
#' @param group_functions Logical. Whether to group function IDs into higher-level
#'   reference categories (via reference maps included in the package). Defaults to \code{TRUE}.
#'
#'
#' @import data.table
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @importFrom phyloseq sample_data tax_table
#'
#' @examples
#' \dontrun{
#' # Example with KO pathway
#' plots <- fun2tax_box(
#'   stratified_path = "ko_metagenome_out_stratified.tsv",
#'   pathway_sdaa_result = sdaa_result,
#'   physeq = physeq_obj,
#'   group = "Treatment",
#'   taxon_level = "Genus",
#'   pathway = "KO"
#' )
#' # Display first plot
#' print(plots[[1]])
#' }
#'
#'
#' @export
fun2tax_box <- function(
    stratified_path,
    pathway_sdaa_result,
    physeq,
    group,
    taxon_level = "Genus",
    pathway     = c("KO","EC","MetaCyc"),
    value_col   = "norm_taxon_function_contrib",
    wrap_width  = 36,
    group_functions = TRUE
){
  # ---- setup ----
  pathway <- match.arg(pathway)
  eps <- 1e-6

  # helpers (unchanged)
  normalize_fun_ids <- function(x, pathway){
    x <- as.character(x)
    if (pathway == "KO") {
      toupper(gsub("^ko:", "", x, ignore.case = TRUE))
    } else if (pathway == "EC") {
      x <- gsub("^ec:", "", x, ignore.case = TRUE)
      toupper(gsub("^EC[\\s_]*", "", x, ignore.case = TRUE))
    } else { # MetaCyc
      toupper(gsub("^metacyc:", "", x, ignore.case = TRUE))
    }
  }
  is_bad_lvl <- function(z){
    z <- as.character(z)
    is.na(z) | !nzchar(z) | z %in% c("Unassigned","Unclassified","Unknown")
  }

  # ---- load stratified table (PICRUSt2) ----
  df <- data.table::fread(stratified_path)
  setDT(df)
  setnames(df, tolower(names(df)))
  # ensure required columns exist (case-insensitive for value_col)
  req <- c("function","taxon","sample")
  if (!all(req %in% names(df))) {
    stop("Stratified file must include columns: function, taxon, sample, and ", value_col)
  }
  if (!value_col %in% names(df)) {
    if (tolower(value_col) %in% names(df)) value_col <- tolower(value_col)
  }
  if (!value_col %in% names(df)) {
    stop("Column '", value_col, "' not found in stratified file.")
  }

  # ---- sample metadata & taxonomy from phyloseq ----
  meta_df <- as.data.frame(phyloseq::sample_data(physeq))
  meta_df$sample <- rownames(meta_df)
  meta <- as.data.table(meta_df)
  setDT(meta)
  if (!group %in% names(meta)) {
    stop("Grouping variable '", group, "' not found in sample_data(physeq).")
  }
  sample_ids <- as.character(meta$sample)

  tax_tab <- as.data.frame(phyloseq::tax_table(physeq))
  tax_tab$ASV <- rownames(tax_tab)
  taxtab <- as.data.table(tax_tab)

  # ---- reference mapping (optional) ----
  load_reference_map <- function(pathway) {
    file_name <- switch(pathway,
                        KO      = "KO_reference.RData",
                        EC      = "EC_reference.RData",
                        MetaCyc = "MetaCyc_reference.RData"
    )
    ref_path <- system.file("extdata", file_name, package = "MicroFun")
    if (!nzchar(ref_path) || !file.exists(ref_path)) {
      stop("Reference file not found in MicroFun/extdata: ", file_name)
    }
    env <- new.env(parent = emptyenv())
    load(ref_path, envir = env)
    ref <- NULL
    for (nm in ls(env)) {
      cand <- get(nm, envir = env)
      if (is.data.frame(cand)) { ref <- cand; break }
    }
    if (is.null(ref)) stop("No data.frame found inside ", file_name)

    if (pathway == "KO") {
      setNames(as.character(ref$PathwayL2),
               as.character(normalize_fun_ids(ref$KO, "KO")))
    } else if (pathway == "EC") {
      setNames(as.character(ref$description),
               as.character(normalize_fun_ids(ref$EC, "EC")))
    } else {
      id_col <- if ("MetaCyc" %in% names(ref)) "MetaCyc" else "ID"
      ids_norm <- normalize_fun_ids(ref[[id_col]], "MetaCyc")
      setNames(as.character(ref$description), as.character(ids_norm))
    }
  }

  id2group <- NULL
  if (isTRUE(group_functions)) {
    id2group <- load_reference_map(pathway)
  }

  # ---- function IDs & grouping labels ----
  fid_vec <- names(pathway_sdaa_result)
  if (is.null(fid_vec) || !length(fid_vec)) {
    stop("pathway_sdaa_result must be a named list keyed by function IDs.")
  }

  grp_vec <- if (is.null(id2group)) {
    rep("All functions", length(fid_vec))
  } else {
    unname(id2group[ normalize_fun_ids(fid_vec, pathway) ])
  }
  grp_vec[is.na(grp_vec) | !nzchar(grp_vec)] <- "Unmapped"

  # ----------------------- FAST PATH STARTS HERE -----------------------
  # 1) Pre-normalize function IDs ONCE for all rows
  df[, fid_norm := normalize_fun_ids(`function`, pathway)]

  # 2) Attach target taxonomic level ONCE via keyed join
  taxtab[, level := .SD[[1]], .SDcols = taxon_level]
  setkey(taxtab, ASV)
  # ensure types
  df[, sample := as.character(sample)]
  # restrict to samples in physeq early
  df <- df[sample %chin% sample_ids]

  # helpful keys
  setkey(df, taxon)           # df$taxon matches taxtab$ASV
  # bring level into df, drop bad levels
  df_annot <- taxtab[df, on = .(ASV = taxon), nomatch = 0L][
    !is_bad_lvl(level)
  ]

  # 3) Materialize 'value' once to avoid get()/SD in groups
  df_annot[, value := .SD[[1]], .SDcols = value_col]

  # 4) Single aggregation to a compact LONG table (fid_norm × level × sample)
  df2 <- df_annot[
    , .(value = sum(value, na.rm = TRUE)),
    by = .(fid_norm, level, sample)
  ]
  # indices for fast lookups
  setindexv(df2, c("fid_norm","level","sample"))
  setindexv(meta, "sample")

  # 5) Normalize the fid_vec once and build groups on normalized IDs
  fid_vec_norm <- normalize_fun_ids(fid_vec, pathway)
  names(fid_vec_norm) <- fid_vec
  groups_norm <- split(fid_vec_norm, grp_vec)

  # 6) Helper: sum across a set of normalized fids → long table
  build_long_for_fids_fast <- function(fids_norm) {
    if (!length(fids_norm)) return(NULL)
    dt <- df2[fid_norm %chin% fids_norm,
              .(abundance = sum(value, na.rm = TRUE)),
              by = .(level, sample)]
    if (!nrow(dt)) return(NULL)
    dt[]
  }

  # ----------------------- PLOTTING (lightweight) ----------------------
  plots <- vector("list", length(groups_norm))
  names(plots) <- names(groups_norm)

  for (glab in names(groups_norm)) {
    fids_n <- groups_norm[[glab]]
    long_sum <- build_long_for_fids_fast(fids_n)
    if (is.null(long_sum) || !nrow(long_sum)) next

    # attach metadata and keep rows that have the grouping variable
    long_sum <- meta[long_sum, on = "sample", nomatch = 0L]
    if (!group %in% names(long_sum)) next
    # long_sum <- long_sum[!is.na(get(group))]
    if (!nrow(long_sum)) next

    # cosmetics
    long_sum[, level_wrapped := stringr::str_wrap(as.character(level), width = wrap_width)]
    long_sum[, abundance_log := log10(pmax(abundance, 0) + eps)]
    long_sum[[group]] <- as.factor(long_sum[[group]])

    p <- ggplot(
      long_sum,
      aes(x = level_wrapped, y = abundance_log, fill = .data[[group]])
    ) +
      geom_boxplot(outlier.shape = 16, outlier.size = 0.8) +
      labs(
        title = glab,
        x = taxon_level,
        y = "Log10 contribution",
        fill = group
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title.position = "plot"
      )

    plots[[glab]] <- p
  }

  Filter(Negate(is.null), plots)
}
