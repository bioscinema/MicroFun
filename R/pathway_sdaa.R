#' Stratified pathway-level DAA workflow for PICRUSt2 outputs
#'
#' Build per-\code{taxon} or per-\code{function} abundance tables from a
#' stratified PICRUSt2 file and run differential abundance analysis (DAA),
#' with optional annotation to KEGG/EC/MetaCyc names. Supports adjusting
#' the modeling design with confounders via backends that accept them
#' (e.g., MaAsLin2).
#'
#' @description Given a PICRUSt2 stratified table (columns \code{sample}, \code{taxon},
#' \code{function}, \code{taxon_function_abun}) and a \code{phyloseq}
#' object, the function:
#' \enumerate{
#'   \item validates overlap of samples and taxa,
#'   \item aggregates ASV-level contributions to a chosen taxonomic rank
#'         (\code{taxon_level}) \emph{or} groups by each \code{function},
#'   \item constructs wide matrices (rows = function IDs; cols = samples),
#'   \item filters tables to ensure at least two groups and sufficient samples,
#'   \item runs \code{pathway_daa()} on each table,
#'   \item annotates results using \code{MicroFun::pathway_annotation()}
#'         or \code{result_annotation()} depending on \code{direction/pathway}.
#' }
#'
#' @param stratified_path Optional character scalar. Path to the stratified PICRUSt2 file
#'   (tab-delimited) containing columns \code{sample}, \code{taxon},
#'   \code{function}, and \code{taxon_function_abun}. For inputs without this
#'   abundance column, \code{norm_taxon_function_contrib} is accepted for
#'   non-KEGG analyses with a compatible DAA method.
#'   Supply exactly one of \code{stratified_path} and \code{stratified_data}.
#' @param stratified_data Optional data frame or data.table containing the same
#'   columns as a contribution file. This accepts the bundled
#'   \code{microfun_demo_KO}, \code{microfun_demo_EC}, \code{microfun_demo_KEGG}
#'   and \code{microfun_demo_MetaCyc} datasets without writing a temporary file.
#'   The supplied object is not modified.
#'
#' @param physeq A \code{phyloseq} object with both \code{tax_table} (ASV rows)
#'   and \code{sample_data}.
#'
#' @param group Character scalar. Column name in \code{sample_data(physeq)} for
#'   the primary grouping variable. Coerced to factor; optionally releveled via
#'   \code{reference}.
#'
#' @param daa_method Character scalar. DAA backend passed to \code{pathway_daa()}
#'   (e.g., \code{'Maaslin2'}, \code{'edgeR'}, \code{'DESeq2'}, \code{'limma voom'},
#'   \code{'metagenomeSeq'}, \code{'Lefser'}).
#'
#' @param p.adjust.method Character scalar. Multiple-testing adjustment method
#'   forwarded to \code{pathway_daa()} (e.g., \code{'none'}, \code{'BH'}).
#'
#' @param direction Character. One of \code{'taxa'} or \code{'function'}.
#'   If \code{'taxa'}, ASV contributions are aggregated to \code{taxon_level}
#'   and analyzed per-taxon table. If \code{'function'}, tables are built per
#'   function (KO/EC/MetaCyc) across taxa.
#'
#' @param pathway Character. Functional namespace of \code{function} IDs:
#'   \code{'KO'}, \code{'EC'}, \code{'MetaCyc'}, or \code{'KEGG'}.
#'   With \code{direction = 'taxa'}, KO input is summarized by KEGG membership,
#'   and KEGG input uses PICRUSt2 pathway contributions directly.
#'
#' @param confounders Optional. Confounding covariates to adjust for in DAA
#'   with MaAsLin2 only. May be:
#'   \itemize{
#'     \item a character vector of column names found in \code{sample_data(physeq)},
#'     \item a numeric/character vector (single confounder),
#'     \item or a \code{data.frame} of covariates (one column per variable).
#'   }
#'   \emph{Note:} This argument is prepared here; use depends on the selected
#'   analysis backend.
#'
#' @param reference Optional character. Reference level for \code{group}; if
#'   present, \code{group} is releveled accordingly.
#'
#' @param min_samples Integer. Require strictly more than \code{min_samples}
#'   samples (after metadata filtering) in a per-table subset to run DAA.
#'   Default: \code{5}.
#'
#' @param verbose Logical. Controls some progress messages. Per-table progress and backend messages may still be printed when FALSE.
#'
#' @param kegg_map Optional PICRUSt2 mapping TSV or named member list passed to
#'   \code{ko2kegg_abundance()} when summarizing KO input.
#'
#' @importFrom data.table dcast fread copy
#'
#' @importFrom reshape2 melt
#'
#' @importFrom phyloseq sample_names tax_table sample_data
#'
#' @importFrom stats aggregate relevel
#'
#' @export
#'
#' @param taxon_level Taxonomic rank in `tax_table(physeq)` used to aggregate ASVs in both orientations; default `"Genus"`.
#'
#' @param maaslin2_transform Transformation used only with `daa_method = "Maaslin2"`: `"AST"` (default), `"LOG"`, `"LOGIT"` or `"NONE"`. MaAsLin2 uses TSS normalization before this transformation. Supplying this argument with another method raises an error.
#'
#' @param threshold Adjusted p-value cutoff passed to annotation. Function-oriented results retain `p_adjust < threshold`. Taxon-oriented KEGG annotation also filters by this cutoff; local KO/EC/MetaCyc annotation does not itself filter p-values.
#'
#' @return A named list of differential-result data frames, keyed by taxon for `direction = "taxa"` or by function for `direction = "function"`. Annotation failures are reported and skipped. The `microfun` attribute records rank, direction, group, namespace, sample IDs, analysis method, reference, adjustment method, contribution column, KEGG mapping and MaAsLin2 transformation. An error is raised if no eligible or retained result remains.
#'
#' @details With `direction = "taxa"` and `pathway = "KO"`, member KO contributions are summed into KEGG pathways before analysis; recorded `function_level` is therefore `"KEGG"`. Direct KEGG input requires `pathway = "KEGG"` and taxon-oriented analysis. These cases require raw `taxon_function_abun`. Other namespaces may use `norm_taxon_function_contrib` with a compatible method. Count-based methods reject fractional input in this workflow. Confounders are supported only with MaAsLin2. Every sample in `physeq` must occur in the contribution file; subset `physeq` explicitly if needed. Results from backends lacking usable p-values or effects may not support downstream annotation or plots.
#'
#' @seealso [pathway_daa()], [tax2fun_box()], [fun2tax_box()]
#'
#' @md
pathway_sdaa <- function(stratified_path = NULL, physeq, taxon_level = "Genus", group, daa_method = "Maaslin2",
  p.adjust.method = "none", threshold = 0.05, direction = "taxa", pathway = "KO", confounders = NULL,
  reference = NULL, min_samples = 5, verbose = TRUE, kegg_map = NULL, maaslin2_transform = "AST",
  stratified_data = NULL) {
  if (daa_method == "Maaslin2") {
    maaslin2_transform <- .microfun_maaslin2_transform(maaslin2_transform)
  } else if (!missing(maaslin2_transform)) {
    stop("maaslin2_transform applies only when daa_method is 'Maaslin2'.")
  }
  stopifnot(direction %in% c("taxa", "function"))
  stopifnot(pathway %in% c("KO", "EC", "MetaCyc", "KEGG"))
  if (pathway == "KEGG" && direction != "taxa")
    stop("Use direction = 'taxa' for KEGG pathway contributions.")
  map_ko <- direction == "taxa" && pathway == "KO"

  # Read stratified contributions and sample metadata.
  df <- .microfun_read_contributions(stratified_path, stratified_data)
  req_cols <- c("sample", "taxon", "function")
  miss <- setdiff(req_cols, names(df))
  if (length(miss))
    stop("Missing required columns in PICRUSt2 file: ", paste(miss, collapse = ", "))
  value_col <- .microfun_contribution_column(names(df), require_abundance = map_ko || pathway == "KEGG")
  values <- df[[value_col]]
  if (!is.numeric(values) || any(!is.finite(values) | values < 0))
    stop("PICRUSt2 contributions must be finite, non-negative numbers.")
  rounding_methods <- c("ALDEx2", "DESeq2", "edgeR", "metagenomeSeq")
  if (value_col == "norm_taxon_function_contrib" &&
    daa_method %in% c(rounding_methods, "limma voom", "LinDA"))
    stop(daa_method, " requires compatible count input; PICRUSt2 proportions or fractional contributions cannot be rounded for this analysis.")

  # phyloseq pieces
  sample_names <- phyloseq::sample_names(physeq)
  if (!length(intersect(as.character(df$sample), sample_names))) {
    stop("No overlapping samples between PICRUSt2 output and phyloseq object.")
  }

  if (!all(sample_names %in% as.character(df$sample)))
    stop("Some physeq samples are completely absent from the stratified file. Subset physeq explicitly.")
  if (map_ko) {
    kegg_map <- .microfun_kegg_members(kegg_map)
    df[["function"]] <- .microfun_ko_ids(df[["function"]])
    df <- df[df[["function"]] %in% unlist(kegg_map, use.names = FALSE), ]
    if (!nrow(df))
      stop("No KO IDs in the stratified file match the KEGG mapping.")
  } else if (pathway == "KEGG") {
    df[["function"]] <- .microfun_kegg_ids(df[["function"]])
    df <- df[grepl("^ko[0-9]{5}$", df[["function"]]), ]
    if (!nrow(df))
      stop("No KEGG pathway IDs found in the stratified file.")
  }
  if (daa_method %in% rounding_methods) {
    values <- df[[value_col]][df$sample %in% sample_names]
    if (any(abs(values - round(values)) > sqrt(.Machine$double.eps)))
      stop(daa_method, " requires integer count input; fractional PICRUSt2 contributions cannot be rounded for this analysis.")
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
      if (length(bad))
        stop("Confounder columns not found in sample_data: ", paste(bad, collapse = ", "))
      conf_df <- mysam[, confounders, drop = FALSE]
    } else if (is.vector(confounders)) {
      conf_df <- data.frame(confounder = confounders)
    } else if (is.data.frame(confounders)) {
      conf_df <- confounders
    } else {
      stop("'confounders' must be NULL, character vector, vector, or data.frame.")
    }
  }

  conf_names <- NULL
  if (!is.null(conf_df)) {
    if (daa_method != "Maaslin2")
      stop("This stratified workflow currently supports confounders only with Maaslin2.")
    if (nrow(conf_df) != nrow(mysam))
      stop("Confounders must have one row per physeq sample.")
    if (is.character(confounders)) {
      conf_names <- confounders
    } else {
      if (identical(rownames(conf_df), as.character(seq_len(nrow(conf_df))))) {
        warning("Unnamed confounders are interpreted in physeq sample order.")
      } else {
        if (!setequal(rownames(conf_df), rownames(mysam)))
          stop("Confounder sample IDs do not match physeq.")
        conf_df <- conf_df[rownames(mysam), , drop = FALSE]
      }
      conf_names <- paste0("MicroFunCovariate", seq_len(ncol(conf_df)))
      mysam[conf_names] <- conf_df
    }
  }
  if (!is.null(reference) && !reference %in% levels(mysam[[group]]))
    stop("reference is not a group level.")

  if (!(taxon_level %in% colnames(taxtab))) {
    stop("The specified 'taxon_level' (", taxon_level, ") is not found in the taxonomy table.")
  }

  # Aggregate contributions by taxon or function.
  if (verbose)
    message("Preparing per-", direction, " tables...")

  make_wide <- function(dt) {
    data.table::dcast(dt, `function` ~ sample, value.var = value_col, fun.aggregate = sum, fill = 0)
  }

  keep_valid_samples <- function(tab, sample_names, id_col = NULL) {
    if (is.null(tab) || NCOL(tab) == 0)
      return(NULL)

    # choose ID column (first by default)
    if (is.null(id_col))
      id_col <- colnames(tab)[1]

    # ensure character and comparable
    sample_names <- as.character(sample_names)

    sample_cols <- intersect(colnames(tab), sample_names)
    cols <- unique(c(id_col, sample_cols))

    if (length(cols) == 1L) {
      # no sample overlap → return NULL to drop this ASV early
      return(NULL)
    }

    if (inherits(tab, "data.table")) {
      # data.table-safe column selection
      tab[, ..cols]
    } else {
      tab[, cols, drop = FALSE]
    }
  }

  # filter taxa (unknown/uncultured) helper
  is_bad_lvl <- function(x) {
    x <- tolower(as.character(x))
    is.na(x) | x == "" | x %in% c("unknown", "uncultured")
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
      make_wide(tmp)
    })

    # Keep only sample columns known to physeq
    per_asv <- lapply(per_asv, keep_valid_samples, sample_names = sample_names)
    per_asv <- Filter(Negate(is.null), per_asv)

    # Merge ASVs to the specified taxon_level (sum by level) Build named list for levels
    if (verbose)
      message("Aggregating ASV tables to level: ", taxon_level)
    result_list_merged <- list()

    for (asv in names(per_asv)) {
      lvl <- map_level[match(asv, map_level$ASV), taxon_level]
      if (length(lvl) == 0 || is_bad_lvl(lvl))
        next
      lvl <- as.character(lvl)[1]

      tab <- per_asv[[asv]]
      if (is.null(tab))
        next
      tab <- as.data.frame(tab, check.names = FALSE)

      # need at least one sample column beyond the id
      if (NCOL(tab) <= 1)
        next
      if (!"function" %in% names(tab))
        next

      # melt + sum + cast ensures aligned union across ASVs within same level
      long <- reshape2::melt(tab, id.vars = "function", variable.name = "sample", value.name = "value")

      # combine into level bucket
      if (is.null(result_list_merged[[lvl]])) {
        result_list_merged[[lvl]] <- long
      } else {
        result_list_merged[[lvl]] <- rbind(result_list_merged[[lvl]], long)
      }
    }

    # sum within level and recast to wide
    for (lvl in names(result_list_merged)) {
      agg <- stats::aggregate(value ~ `function` + sample, data = result_list_merged[[lvl]], sum)
      wide <- reshape2::dcast(agg, `function` ~ sample, value.var = "value", fill = 0)
      result_list_merged[[lvl]] <- wide
    }

    # Filter by metadata availability and sample counts; store
    for (lvl in names(result_list_merged)) {
      tbl <- result_list_merged[[lvl]]
      # intersect with mysam samples (drop those missing group)
      sample_ids <- setdiff(colnames(tbl), "function")
      sample_ids <- rownames(mysam)
      for (sid in setdiff(sample_ids, colnames(tbl))) tbl[[sid]] <- 0
      g <- mysam[[group]][match(sample_ids, rownames(mysam))]
      ok <- !is.na(g)
      sample_ids <- sample_ids[ok]
      g <- g[ok]

      if (length(unique(g)) < 2) {
        if (verbose)
          message("Skipping ", lvl, " (only one group present).")
        next
      }
      if (length(sample_ids) <= min_samples) {
        if (verbose)
          message("Skipping ", lvl, " (", length(sample_ids), " samples; need >", min_samples,
          ").")
        next
      }
      result_list_subset[[lvl]] <- tbl[, c("function", sample_ids), drop = FALSE]
    }

  } else {
    # Split by function (KO/EC/MetaCyc)
    fun_list <- split(df, df$`function`)
    for (fname in names(fun_list)) {
      dt <- fun_list[[fname]]
      # remove unclassified ASV
      dt <- dt[dt$taxon != "unclassified", ]
      # map ASV -> level
      dt$level <- taxtab[match(dt$taxon, taxtab$ASV), taxon_level]
      dt <- dt[!is_bad_lvl(dt$level), ]

      if (!nrow(dt))
        next
      # aggregate by level x sample
      dt_agg <- dt[, .(value = sum(get(value_col), na.rm = TRUE)),
        by = .(taxon = level, sample)]
      wide <- data.table::dcast(dt_agg, taxon ~ sample, value.var = "value",
        fill = 0)

      # keep only samples with metadata
      sample_ids <- intersect(colnames(wide), sample_names)
      present <- rownames(mysam)
      for (sid in setdiff(present, colnames(wide))) wide[, (sid) := 0]
      g <- mysam[[group]][match(present, rownames(mysam))]
      ok <- !is.na(g)
      present <- present[ok]
      g <- g[ok]

      if (length(unique(g)) < 2)
        next
      if (length(present) <= min_samples)
        next

      result_list_subset[[fname]] <- wide[, c("taxon", present), with = FALSE]
    }
  }

  if (!length(result_list_subset)) {
    stop("No eligible ", direction, " tables after filtering by groups and sample counts.")
  }
  run_daa <- function(abundance, metadata) {
    args <- list(abundance = abundance, metadata = metadata, group = group, daa_method = daa_method,
      p.adjust = p.adjust.method, confounder = conf_names, reference = reference)
    if (daa_method == "Maaslin2")
      args$maaslin2_transform <- maaslin2_transform
    do.call(pathway_daa, args)
  }

  # Run differential analysis and annotate each result.
  if (verbose)
    message("Running DAA on ", length(result_list_subset), " ", direction, " tables...")

  # Initialize result list
  all_daa_results <- list()

  # Iterate through taxon-level KO tables with a progress bar
  total <- length(result_list_subset)
  for (i in seq_along(result_list_subset)) {
    name <- names(result_list_subset)[i]
    cat(sprintf("\n[%d/%d] Processing %s: %s\n", i, total, ifelse(direction == "taxa", "taxon", "function"),
      name))

    cat("Running differential analysis: [", sep = "")
    for (j in 1:30) {
      cat("=")
      flush.console()
      Sys.sleep(0.01)
    }
    cat("]\n")

    df1 <- data.frame(result_list_subset[[name]], check.names = FALSE)

    if (direction == "taxa") {
      if (map_ko) {
        ko_to_kegg <- TRUE

        # Sum predicted contributions for the member KOs of each pathway.
        if (!any(rowSums(df1[, -1, drop = FALSE]) > 0))
          next
        ko_abun <- ko2kegg_abundance(data = df1, kegg_map = kegg_map)

      } else if (pathway == "EC") {
        # Remove rows with NA or empty function values
        df1 <- df1[!is.na(df1$`function`) & df1$`function` != "", ]

        rownames(df1) <- paste0(pathway, df1$`function`)
        ko_abun <- df1[, setdiff(colnames(df1), "function"), drop = FALSE]
        ko_to_kegg <- FALSE

      } else if (pathway %in% c("MetaCyc", "KEGG")) {
        # Remove rows with NA or empty function values
        df1 <- df1[!is.na(df1$`function`) & df1$`function` != "", ]

        rownames(df1) <- df1$`function`
        ko_abun <- df1[, setdiff(colnames(df1), "function"), drop = FALSE]
        ko_to_kegg <- pathway == "KEGG"
      }

      # Prepare metadata
      mysam_subset <- mysam[rownames(mysam) %in% colnames(ko_abun), ]
      mysam_subset$sample <- rownames(mysam_subset)
      mysam_subset[[group]] <- as.factor(mysam_subset[[group]])

      # DAA
      result <- run_daa(ko_abun, mysam_subset)
      if ("feature" %in% colnames(result)) {
        result <- result[!is.na(result$feature) & result$feature != "", ]
      }

      # Annotate and store result
      annotated <- tryCatch({
        MicroFun::pathway_annotation(pathway = pathway, daa_results_df = result, ko_to_kegg = ko_to_kegg,
          threshold = threshold)
      }, error = function(e) {
        message("  -> Skipping ", name, ": ", e$message)
        NULL
      })

    } else if (direction == "function") {
      rownames(df1) <- df1[, "taxon"]
      df1 <- df1[, -which(colnames(df1) == "taxon")]

      # Prepare metadata
      mysam_subset <- mysam[rownames(mysam) %in% colnames(df1), ]
      mysam_subset$sample <- rownames(mysam_subset)
      mysam_subset[[group]] <- as.factor(mysam_subset[[group]])

      # DAA directly on function-level abundance
      result <- run_daa(df1, mysam_subset)
      if ("feature" %in% colnames(result)) {
        result <- result[!is.na(result$feature) & result$feature != "", ]
      }

      annotated <- tryCatch({
        res <- result_annotation(result, threshold = threshold)
        if (is.null(res) || nrow(res) == 0) {
          message("  -> Skipping ", name, ": no significant features found")
          res <- NULL
        }
        res
      }, error = function(e) {
        message("  -> Skipping ", name, ": ", e$message)
        NULL
      })

    }

    # Store the analysis result.
    if (!is.null(annotated)) {
      all_daa_results[[name]] <- annotated
    }
  }
  if (length(all_daa_results) == 0) {
    stop(sprintf("No statistically significant biomarkers found for any %s.\nTry a less stringent threshold or review your input data.",
      ifelse(direction == "taxa", "taxon", "function")))
  }
  attr(all_daa_results, "microfun") <- list(taxon_level = taxon_level, direction = direction, group = group,
    pathway = pathway, function_level = if (map_ko) "KEGG" else pathway,
    samples = rownames(mysam), daa_method = daa_method, reference = reference, p.adjust.method = p.adjust.method,
    value_col = value_col, kegg_map = if (map_ko) kegg_map,
    maaslin2_transform = if (daa_method == "Maaslin2") maaslin2_transform)
  return(all_daa_results)
}

# Keep the reader separate so file and in-memory inputs share all downstream
# validation and analysis. Copy tables before any data.table mutation.
.microfun_read_contributions <- function(stratified_path = NULL, stratified_data = NULL) {
  if (is.null(stratified_path) == is.null(stratified_data)) {
    stop("Supply exactly one of stratified_path or stratified_data.")
  }
  if (!is.null(stratified_data)) {
    if (!is.data.frame(stratified_data)) {
      stop("stratified_data must be a data frame or data.table.")
    }
    return(data.table::copy(data.table::as.data.table(stratified_data)))
  }
  if (!is.character(stratified_path) || length(stratified_path) != 1L ||
      is.na(stratified_path) || !nzchar(stratified_path) || !file.exists(stratified_path)) {
    stop("stratified_path must name an existing contribution file.")
  }
  if (grepl("\\.gz$", stratified_path, ignore.case = TRUE) &&
      !requireNamespace("R.utils", quietly = TRUE)) {
    con <- gzfile(stratified_path, "rt")
    on.exit(close(con))
    return(data.table::as.data.table(utils::read.delim(con, check.names = FALSE)))
  }
  data.table::fread(stratified_path, sep = "\t")
}
