#' Stratified pathway-level DAA from PICRUSt2 outputs
#'
#' Build per-\code{taxon} or per-\code{function} abundance tables from a
#' stratified PICRUSt2 file and run differential abundance analysis (DAA)
#' at the pathway/function level, optionally adjusting for confounders.
#'
#' @description
#' Given a PICRUSt2 stratified output (e.g., \code{pred_metagenome_contrib.tsv})
#' and a \code{phyloseq} object, this function:
#' \enumerate{
#'   \item checks sample and taxa overlap,
#'   \item aggregates ASV-level contributions to a chosen taxonomy level
#'         (e.g., \code{"Genus"}) \emph{or} groups by function,
#'   \item constructs wide matrices (rows = KO/EC/MetaCyc IDs or taxa, cols = samples),
#'   \item filters tables with insufficient group representation,
#'   \item runs DAA via \code{pathway_daa()} (supports edgeR/DESeq2/limma-voom/Maaslin2/etc.),
#'   \item annotates significant results via \code{MicroFun::pathway_annotation()} or
#'         \code{result_annotation()}.
#' }
#'
#' @param stratified_path \strong{character}. Path to PICRUSt2 stratified output
#'   (tab-delimited) containing at least the columns:
#'   \code{sample}, \code{taxon}, \code{function},
#'   \code{norm_taxon_function_contrib}.
#' @param physeq \strong{phyloseq} object. Must contain both \code{tax_table}
#'   (with an \code{ASV} rownames key) and \code{sample_data}.
#' @param taxon_level \strong{character}. Taxonomic rank to aggregate to when
#'   \code{direction = "taxa"} (e.g., \code{"Genus"}, \code{"Family"}). Must be a
#'   column in \code{phyloseq::tax_table(physeq)}.
#' @param group \strong{character}. Column name in \code{sample_data} indicating
#'   the primary grouping variable for DAA (will be coerced to a factor).
#' @param daa_method \strong{character}. DAA backend name passed to
#'   \code{pathway_daa()} (e.g., \code{"Maaslin2"}, \code{"edgeR"},
#'   \code{"DESeq2"}, \code{"limma voom"}, \code{"metagenomeSeq"}, \code{"Lefser"}).
#' @param p.adjust.method \strong{character}. Multiple-testing correction method
#'   forwarded to \code{pathway_daa()} (e.g., \code{"none"}, \code{"BH"}).
#' @param threshold \strong{numeric}. Significance cutoff used during annotation
#'   (typically FDR/q-value threshold).
#' @param direction \strong{character}. One of \code{"taxa"} or \code{"function"}.
#'   If \code{"taxa"}, the function aggregates ASV contributions to \code{taxon_level}
#'   and analyzes per-taxon tables. If \code{"function"}, it groups by function
#'   (KO/EC/MetaCyc) and analyzes per-function tables.
#' @param pathway \strong{character}. Functional namespace of the \code{function}
#'   column: one of \code{"KO"}, \code{"EC"}, or \code{"MetaCyc"}. Controls
#'   downstream annotation and (for KO) \code{ko_to_kegg} mapping.
#' @param confounders Optional covariate(s) to adjust for in DAA backends that
#'   support them (e.g., edgeR/limma/Maaslin2). Either a vector (single covariate)
#'   or a \code{data.frame}/tibble with one column per covariate. If a character
#'   vector, the named columns are pulled from \code{sample_data(physeq)}. Must
#'   align with samples in the abundance matrices. Characters are converted to
#'   factors by the modeling backend.
#' @param reference \strong{character} or \code{NULL}. Optional reference level
#'   for \code{group}; if supplied and present, \code{group} is releveled accordingly.
#' @param min_samples \strong{integer}. Require strictly more than this number of
#'   samples in the per-table subset (after filtering to samples with non-missing
#'   group) to run DAA. Default: \code{5}.
#' @param verbose \strong{logical}. If \code{TRUE}, print progress messages.
#'
#' @examples
#' \dontrun{
#' # Example inputs (toy):
#' stratified_path <- "pred_metagenome_contrib.tsv"
#' group <- "Group"
#'
#' # 1) Per-taxon KO analysis at Genus level using Maaslin2
#' res_taxa <- pathway_sdaa(
#'   stratified_path = stratified_path,
#'   physeq          = physeq_obj,
#'   taxon_level     = "Genus",
#'   group           = group,
#'   daa_method      = "Maaslin2",
#'   p.adjust.method = "BH",
#'   threshold       = 0.05,
#'   direction       = "taxa",
#'   pathway         = "KO",
#'   confounders     = c("age","sex"),
#'   reference       = "Control"
#' )
#'
#' # 2) Per-function MetaCyc analysis without confounders
#' res_fun <- pathway_sdaa(
#'   stratified_path = stratified_path,
#'   physeq          = physeq_obj,
#'   taxon_level     = "Genus",
#'   group           = group,
#'   daa_method      = "edgeR",
#'   p.adjust.method = "BH",
#'   threshold       = 0.1,
#'   direction       = "function",
#'   pathway         = "MetaCyc"
#' )
#' }
#'
#' @export
#'
pathway_sdaa <- function(
    stratified_path,
    physeq,
    taxon_level = "Genus",
    group,
    daa_method = "Maaslin2",
    p.adjust.method = "none",
    threshold = 0.05,
    direction = "taxa",          # "taxa" or "function"
    pathway   = "KO",            # "KO", "EC", or "MetaCyc"
    confounders = NULL,          # NULL, character vector of metadata columns, vector, or data.frame
    reference   = NULL,          # optional reference level for `group`
    min_samples = 5,             # require > min_samples per test table
    verbose     = TRUE
) {
  stopifnot(direction %in% c("taxa","function"))
  stopifnot(pathway %in% c("KO","EC","MetaCyc"))

  # --- 1) Load & validate inputs ------------------------------------------------
  df <- data.table::fread(stratified_path, sep = "\t")
  req_cols <- c("sample","taxon","function","norm_taxon_function_contrib")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing required columns in PICRUSt2 file: ", paste(miss, collapse = ", "))

  # phyloseq pieces
  sample_names <- phyloseq::sample_names(physeq)
  if (!length(intersect(as.character(df$sample), sample_names))) {
    stop("No overlapping samples between PICRUSt2 output and phyloseq object.")
  }

  taxtab <- as.data.frame(phyloseq::tax_table(physeq))
  taxtab$ASV <- rownames(taxtab)

  # sample_data as data.frame with rownames as sample ids
  mysam <- as.data.frame(phyloseq::sample_data(physeq))
  if (!(group %in% colnames(mysam))) {
    stop("Group column '", group, "' is not present in sample_data.")
  }

  # enforce factor for group, optionally relevel by reference
  mysam[[group]] <- as.factor(mysam[[group]])
  if (!is.null(reference) && reference %in% levels(mysam[[group]])) {
    mysam[[group]] <- stats::relevel(mysam[[group]], ref = reference)
  }

  # prepare confounders (if character, pull from mysam)
  conf_df <- NULL
  if (!is.null(confounders)) {
    if (is.character(confounders)) {
      bad <- setdiff(confounders, colnames(mysam))
      if (length(bad)) stop("Confounder columns not found in sample_data: ", paste(bad, collapse = ", "))
      conf_df <- mysam[, confounders, drop = FALSE]
    } else if (is.vector(confounders)) {
      conf_df <- data.frame(confounder = confounders)
    } else if (is.data.frame(confounders)) {
      conf_df <- confounders
    } else {
      stop("'confounders' must be NULL, character vector, vector, or data.frame.")
    }
  }

  if (!(taxon_level %in% colnames(taxtab))) {
    stop("The specified 'taxon_level' (", taxon_level, ") is not found in the taxonomy table.")
  }

  # --- 2) Build KO/EC/MetaCyc tables per taxon or per function ------------------
  if (verbose) message("Preparing per-", direction, " tables...")

  make_wide <- function(dt, value_col = "norm_taxon_function_contrib") {
    # dt: must have a 'row_id' column for rows and 'sample' for columns
    data.table::dcast(dt, row_id ~ sample, value.var = value_col, fill = 0)
  }

  keep_valid_samples <- function(tab, sample_names) {
    sample_cols <- intersect(colnames(tab), sample_names)
    tab[, c(colnames(tab)[1], sample_cols), drop = FALSE]
  }

  # filter taxa (unknown/uncultured) helper
  is_bad_lvl <- function(x) {
    x <- tolower(as.character(x))
    is.na(x) | x == "" | x %in% c("unknown","uncultured")
  }

  result_list_subset <- list()

  if (direction == "taxa") {
    # Split by ASV/taxon
    taxon_list <- split(df, df$taxon)

    # ASV -> taxon_level mapping
    map_level <- taxtab[, c("ASV", taxon_level)]
    # Per ASV build wide: rows=function (KO/EC/MetaCyc id), cols=samples
    per_asv <- lapply(taxon_list, function(dt) {
      tmp <- data.table::copy(dt)
      tmp[, row_id := `function`]
      make_wide(tmp)
    })

    # Keep only sample columns known to physeq
    per_asv <- lapply(per_asv, keep_valid_samples, sample_names = sample_names)

    # Merge ASVs to the specified taxon_level (sum by level)
    # Build named list for levels
    if (verbose) message("Aggregating ASV tables to level: ", taxon_level)
    result_list_merged <- list()

    for (asv in names(per_asv)) {
      lvl <- map_level[match(asv, map_level$ASV), taxon_level]
      if (length(lvl) == 0 || is_bad_lvl(lvl)) next

      tab <- per_asv[[asv]]
      if (!ncol(tab) > 1) next
      # melt + sum + cast ensures aligned union across ASVs within same level
      long <- reshape2::melt(tab, id.vars = colnames(tab)[1], variable.name = "sample", value.name = "value")
      # combine into level bucket
      if (is.null(result_list_merged[[lvl]])) {
        result_list_merged[[lvl]] <- long
      } else {
        result_list_merged[[lvl]] <- rbind(result_list_merged[[lvl]], long)
      }
    }

    # sum within level and recast to wide
    for (lvl in names(result_list_merged)) {
      agg <- stats::aggregate(value ~ row_id + sample, data = result_list_merged[[lvl]], sum)
      wide <- reshape2::dcast(agg, row_id ~ sample, value.var = "value", fill = 0)
      result_list_merged[[lvl]] <- wide
    }

    # Filter by metadata availability and sample counts; store
    for (lvl in names(result_list_merged)) {
      tbl <- result_list_merged[[lvl]]
      # intersect with mysam samples (drop those missing group)
      sample_ids <- setdiff(colnames(tbl), "row_id")
      sample_ids <- intersect(sample_ids, rownames(mysam))
      g <- mysam[sample_ids, group]
      ok <- !is.na(g)
      sample_ids <- sample_ids[ok]
      g <- g[ok]

      if (length(unique(g)) < 2) {
        if (verbose) message("Skipping ", lvl, " (only one group present).")
        next
      }
      if (length(sample_ids) <= min_samples) {
        if (verbose) message("Skipping ", lvl, " (", length(sample_ids), " samples; need >", min_samples, ").")
        next
      }
      result_list_subset[[lvl]] <- tbl[, c("row_id", sample_ids), drop = FALSE]
    }

  } else { # direction == "function"
    # Split by function (KO/EC/MetaCyc)
    fun_list <- split(df, df$`function`)
    for (fname in names(fun_list)) {
      dt <- fun_list[[fname]]
      # remove unclassified ASV
      dt <- dt[dt$taxon != "unclassified", ]
      # map ASV -> level
      dt$level <- taxtab[match(dt$taxon, taxtab$ASV), taxon_level]
      dt <- dt[!is_bad_lvl(dt$level), ]

      if (!nrow(dt)) next
      # aggregate by level x sample
      dt_agg <- dt[, .(norm_taxon_function_contrib = sum(norm_taxon_function_contrib, na.rm = TRUE)),
                   by = .(row_id = level, sample)]
      wide <- data.table::dcast(dt_agg, row_id ~ sample, value.var = "norm_taxon_function_contrib", fill = 0)

      # keep only samples with metadata
      sample_ids <- intersect(colnames(wide), sample_names)
      present <- intersect(sample_ids, rownames(mysam))
      g <- mysam[present, group]
      ok <- !is.na(g)
      present <- present[ok]; g <- g[ok]

      if (length(unique(g)) < 2) next
      if (length(present) <= min_samples) next

      result_list_subset[[fname]] <- wide[, c("row_id", present), drop = FALSE]
    }
  }

  if (!length(result_list_subset)) {
    stop("No eligible ", direction, " tables after filtering by groups and sample counts.")
  }

  # --- 3) DAA per table + annotation -------------------------------------------
  if (verbose) message("Running DAA on ", length(result_list_subset), " ", direction, " tables...")

  all_daa_results <- list()

  for (nm in names(result_list_subset)) {
    if (verbose) message("Processing ", direction, ": ", nm)

    tbl <- result_list_subset[[nm]]
    # Prepare abundance & metadata
    rownames(tbl) <- tbl$row_id
    tbl <- tbl[, -1, drop = FALSE]

    mysam_subset <- mysam[rownames(mysam) %in% colnames(tbl), , drop = FALSE]
    mysam_subset$sample <- rownames(mysam_subset)
    mysam_subset[[group]] <- droplevels(as.factor(mysam_subset[[group]]))

    # DAA
    # Note: pathway_daa should internally route to edgeR/DESeq2/etc.; we pass confounders & reference when supported
    result <- pathway_daa(
      abundance = tbl,
      metadata  = mysam_subset,
      group     = group,
      daa_method= daa_method,
      p.adjust  = p.adjust.method,
      confounder= conf_df,
      reference = reference
    )

    if (!is.null(result) && "feature" %in% colnames(result)) {
      result <- result[!is.na(result$feature) & result$feature != "", , drop = FALSE]
    }

    # annotation
    annotated <- NULL
    if (direction == "taxa") {
      ko_to_kegg <- FALSE
      if (pathway == "KO") {
        # strip "ko:" prefix if present before ko2kegg_abundance inside DAA step
        # (here only for annotation mapping flags)
        ko_to_kegg <- TRUE
      }
      annotated <- tryCatch(
        MicroFun::pathway_annotation(
          pathway        = pathway,
          daa_results_df = result,
          ko_to_kegg     = ko_to_kegg,
          threshold      = threshold
        ),
        error = function(e) {
          if (verbose) message("  -> Skipping ", nm, ": ", e$message)
          NULL
        }
      )
    } else {
      annotated <- tryCatch(
        {
          res <- result_annotation(result, threshold = threshold)
          if (is.null(res) || !nrow(res)) NULL else res
        },
        error = function(e) {
          if (verbose) message("  -> Skipping ", nm, ": ", e$message)
          NULL
        }
      )
    }

    if (!is.null(annotated) && nrow(annotated)) {
      all_daa_results[[nm]] <- annotated
    }
  }

  if (!length(all_daa_results)) {
    stop(sprintf(
      "No statistically significant biomarkers found for any %s.\nTry a less stringent threshold or review your input data.",
      ifelse(direction == "taxa","taxon","function")
    ))
  }

  all_daa_results
}
