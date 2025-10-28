#' Stratified pathway-level DAA workflow for PICRUSt2 outputs
#'
#' Build per-\code{taxon} or per-\code{function} abundance tables from a
#' stratified PICRUSt2 file and run differential abundance analysis (DAA),
#' with optional annotation to KEGG/EC/MetaCyc names. Supports adjusting
#' the modeling design with confounders via backends that accept them
#' (e.g., MaAsLin2).
#'
#' @description
#' Given a PICRUSt2 stratified table (columns \code{sample}, \code{taxon},
#' \code{function}, \code{norm_taxon_function_contrib}) and a \code{phyloseq}
#' object, the function:
#' \enumerate{
#'   \item validates overlap of samples and taxa,
#'   \item aggregates ASV-level contributions to a chosen taxonomic rank
#'         (\code{taxon_level}) \emph{or} groups by each \code{function},
#'   \item constructs wide matrices (rows = function IDs; cols = samples),
#'   \item filters tables to ensure at least two groups and sufficient samples,
#'   \item runs \code{pathway_daa()} on each table,
#'   \item annotates significant results using \code{MicroFun::pathway_annotation()}
#'         or \code{result_annotation()} depending on \code{direction/pathway}.
#' }
#'
#' @param stratified_path Character scalar. Path to the stratified PICRUSt2 file
#'   (tab-delimited) containing columns \code{sample}, \code{taxon},
#'   \code{function}, \code{norm_taxon_function_contrib}.
#' @param physeq A \code{phyloseq} object with both \code{tax_table} (ASV rows)
#'   and \code{sample_data}.
#' @param taxon_level Character scalar. Taxonomic rank used to aggregate ASVs
#'   when \code{direction = "taxa"} (e.g., \code{"Genus"}). Must exist in
#'   \code{phyloseq::tax_table(physeq)}.
#' @param group Character scalar. Column name in \code{sample_data(physeq)} for
#'   the primary grouping variable. Coerced to factor; optionally releveled via
#'   \code{reference}.
#' @param daa_method Character scalar. DAA backend passed to \code{pathway_daa()}
#'   (e.g., \code{"Maaslin2"}, \code{"edgeR"}, \code{"DESeq2"}, \code{"limma voom"},
#'   \code{"metagenomeSeq"}, \code{"Lefser"}).
#' @param p.adjust.method Character scalar. Multiple-testing adjustment method
#'   forwarded to \code{pathway_daa()} (e.g., \code{"none"}, \code{"BH"}).
#' @param threshold Numeric. Significance cutoff used during result annotation
#'   (typically FDR/q-value threshold). Default: \code{0.05}.
#' @param direction Character. One of \code{"taxa"} or \code{"function"}.
#'   If \code{"taxa"}, ASV contributions are aggregated to \code{taxon_level}
#'   and analyzed per-taxon table. If \code{"function"}, tables are built per
#'   function (KO/EC/MetaCyc) across taxa.
#' @param pathway Character. Functional namespace of \code{function} IDs:
#'   \code{"KO"}, \code{"EC"}, or \code{"MetaCyc"}. Controls annotation and
#'   KO→KEGG pathway summarization.
#' @param confounders Optional. Confounding covariates to adjust for in DAA
#'   backends that support them. May be:
#'   \itemize{
#'     \item a character vector of column names found in \code{sample_data(physeq)},
#'     \item a numeric/character vector (single confounder),
#'     \item or a \code{data.frame} of covariates (one column per variable).
#'   }
#'   \emph{Note:} This argument is prepared here; actual use depends on your
#'   \code{pathway_daa()} implementation (ensure you forward it).
#' @param reference Optional character. Reference level for \code{group}; if
#'   present, \code{group} is releveled accordingly.
#' @param min_samples Integer. Require strictly more than \code{min_samples}
#'   samples (after metadata filtering) in a per-table subset to run DAA.
#'   Default: \code{5}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#'
#'
#' @examples
#' \dontrun{
#' stratified_path <- "pred_metagenome_contrib.tsv"
#' # Minimal inputs
#' out <- pathway_sdaa(
#'   stratified_path = stratified_path,
#'   physeq          = physeq_obj,
#'   taxon_level     = "Genus",
#'   group           = "Group",
#'   daa_method      = "Maaslin2",
#'   p.adjust.method = "BH",
#'   threshold       = 0.05,
#'   direction       = "taxa",
#'   pathway         = "KO",
#'   confounders     = c("Age","Sex"),
#'   reference       = "Control",
#'   min_samples     = 5
#' )
#'
#' # Function-centric EC analysis:
#' out_fun <- pathway_sdaa(
#'   stratified_path = stratified_path,
#'   physeq          = physeq_obj,
#'   taxon_level     = "Genus",
#'   group           = "Group",
#'   daa_method      = "edgeR",
#'   direction       = "function",
#'   pathway         = "EC"
#' )
#' }
#'
#'
#' @importFrom data.table dcast fread copy
#' @importFrom reshape2 melt
#' @importFrom phyloseq sample_names tax_table sample_data
#' @importFrom stats aggregate relevel
#' @export
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
    data.table::dcast(dt, `function` ~ sample, value.var = value_col, fill = 0)
  }

  keep_valid_samples <- function(tab, sample_names, id_col = NULL) {
    if (is.null(tab) || NCOL(tab) == 0) return(NULL)

    # choose ID column (first by default)
    if (is.null(id_col)) id_col <- colnames(tab)[1]

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
      make_wide(tmp)
    })

    # Keep only sample columns known to physeq
    per_asv <- lapply(per_asv, keep_valid_samples, sample_names = sample_names)
    per_asv <- Filter(Negate(is.null), per_asv)

    # Merge ASVs to the specified taxon_level (sum by level)
    # Build named list for levels
    if (verbose) message("Aggregating ASV tables to level: ", taxon_level)
    result_list_merged <- list()

    for (asv in names(per_asv)) {
      lvl <- map_level[match(asv, map_level$ASV), taxon_level]
      if (length(lvl) == 0 || is_bad_lvl(lvl)) next
      lvl <- as.character(lvl)[1]

      tab <- per_asv[[asv]]
      if (is.null(tab)) next
      tab <- as.data.frame(tab, check.names = FALSE)

      # need at least one sample column beyond the id
      if (NCOL(tab) <= 1) next
      if (!"function" %in% names(tab)) next

      # melt + sum + cast ensures aligned union across ASVs within same level
      long <- reshape2::melt(
        tab,
        id.vars = "function",
        variable.name = "sample",
        value.name = "value"
      )

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
      sample_ids <- intersect(sample_ids, rownames(mysam))
      g <- mysam[sample_ids, group]
      ok <- !is.na(g)
      sample_ids <- sample_ids[ok]
      g <- g[ok]

      if (length(unique(g[[group]])) < 2) {
        if (verbose) message("Skipping ", lvl, " (only one group present).")
        next
      }
      if (length(sample_ids) <= min_samples) {
        if (verbose) message("Skipping ", lvl, " (", length(sample_ids), " samples; need >", min_samples, ").")
        next
      }
      result_list_subset[[lvl]] <- tbl[, c("function", sample_ids), drop = FALSE]
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
                   by = .(taxon = level, sample)]
      wide <- data.table::dcast(dt_agg, taxon ~ sample, value.var = "norm_taxon_function_contrib", fill = 0)

      # keep only samples with metadata
      sample_ids <- intersect(colnames(wide), sample_names)
      present <- intersect(sample_ids, rownames(mysam))
      g <- mysam[present, group]
      ok <- !is.na(g)
      present <- present[ok]; g <- g[ok]

      if (length(unique(g[[group]])) < 2) next
      if (length(present) <= min_samples) next

      result_list_subset[[fname]] <- wide[, c("taxon", present), with = FALSE]
    }
  }

  if (!length(result_list_subset)) {
    stop("No eligible ", direction, " tables after filtering by groups and sample counts.")
  }

  # --- 3) DAA per table + annotation -------------------------------------------
  if (verbose) message("Running DAA on ", length(result_list_subset), " ", direction, " tables...")

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
      if (pathway == "KO") {
        # Remove rows with NA or empty function values
        df1 <- df1[!is.na(df1$`function`) & df1$`function` != "", ]

        # Remove "ko:" prefix from function names
        df1$`function` <- gsub("^ko:", "", df1$`function`)
        ko_to_kegg <- TRUE

        # Convert from KO-level to KEGG pathway abundance
        ko_abun <- ko2kegg_abundance(data = df1)

      } else if (pathway == "EC") {
        # Remove rows with NA or empty function values
        df1 <- df1[!is.na(df1$`function`) & df1$`function` != "", ]

        rownames(df1) <- paste0(pathway, df1$`function`)
        ko_abun <- df1[, setdiff(colnames(df1), "function"), drop = FALSE]
        ko_to_kegg <- FALSE

      } else if (pathway == "MetaCyc") {
        # Remove rows with NA or empty function values
        df1 <- df1[!is.na(df1$`function`) & df1$`function` != "", ]

        rownames(df1) <- df1$`function`
        ko_abun <- df1[, setdiff(colnames(df1), "function"), drop = FALSE]
        ko_to_kegg <- FALSE
      }



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
      if ("feature" %in% colnames(result)) {
        result <- result[!is.na(result$feature) & result$feature != "", ]
      }
      # result is a data.frame/matrix with rownames like "EC:1.1.1.1", "ec_2.7.1.1", "EC1.2.3.4"
      # result$feature <- sub("(?i)^EC[:_\\-]*", "", result$feature, perl = TRUE)


      # Annotate and store result
      annotated <- tryCatch(
        {
          MicroFun::pathway_annotation(
            pathway        = pathway,
            daa_results_df = result,
            ko_to_kegg     = ko_to_kegg,
            threshold      = threshold,

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
      rownames(df1) <- df1[,"taxon"]
      df1 <- df1[, -which(colnames(df1) == "taxon")]

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
      if ("feature" %in% colnames(result)) {
        result <- result[!is.na(result$feature) & result$feature != "", ]
      }

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




