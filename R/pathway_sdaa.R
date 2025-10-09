#' Differential pathway analysis from PICRUSt2 stratified output (by taxon or by function)
#'
#' @description
#' Runs a per-*taxon* or per-*function* differential abundance analysis (DAA) on
#' stratified PICRUSt2 output matched to a `phyloseq` object. Optionally aggregates
#' ASVs to a specified taxonomic rank, converts KO to KEGG pathway abundance when
#' requested, and annotates significant results.
#'
#' @param stratified_path character. Path to the PICRUSt2 *stratified* table
#'   (tab-delimited). The file must contain at least the columns:
#'   `sample`, `taxon`, `function`, and either
#'   `taxon_function_abun` (for `direction = "taxa"`) or
#'   `norm_taxon_function_contrib` (for `direction = "function"`).
#' @param physeq A `phyloseq` object containing `sample_data()` and `tax_table()`
#'   with taxa names matching the PICRUSt2 `taxon` field.
#' @param taxon_level character. Taxonomic rank to merge ASVs to (e.g., `"Genus"`).
#'   If `NULL`, no taxonomic aggregation is performed. Must be a column of
#'   `tax_table(physeq)`. Default `"Genus"`.
#' @param group character (scalar). Column name in `sample_data(physeq)` defining
#'   the groups/phenotypes used in the DAA design.
#' @param daa_method character. Method passed to `pathway_daa()` (e.g., `"Maaslin2"`).
#'   Default `"Maaslin2"`.
#' @param p.adjust.method character. Multiple-testing correction method passed to
#'   `pathway_daa()` (e.g., `"BH"`, `"none"`). Default `"none"`.
#' @param threshold numeric. Significance threshold used by the annotation step
#'   (`result_annotation()` / `MicroFun::pathway_annotation()`). Default `0.05`.
#' @param direction character. Workflow direction: `"taxa"` (split by taxon/ASV,
#'   build function × sample tables per taxon) or `"function"` (split by function
#'   and aggregate contributions by the chosen taxon level). Default `"taxa"`.
#' @param pathway character. Pathway/function namespace:
#'   `"KO"`, `"EC"`, or `"MetaCyc"`.
#'   * `"KO"`: strips `"ko:"` prefix and converts KO → KEGG pathway abundance via
#'   `ko2kegg_abundance()` (sets `ko_to_kegg = TRUE`);
#'   * `"EC"` / `"MetaCyc"`: uses the `function` rows as-is (no KO→KEGG conversion;
#'   sets `ko_to_kegg = FALSE`).
#'
#' @details
#' **Pipeline (high level):**
#' 1. Read stratified table (`data.table::fread`) and verify overlap with
#'    `sample_names(physeq)` and `taxa_names(physeq)`.
#' 2. If `direction = "taxa"`:
#'    - Split by `taxon`; reshape to wide (rows = `function`, cols = samples) using
#'      `taxon_function_abun`.
#'    - Optionally merge ASVs to `taxon_level` using `tax_table(physeq)` and sum.
#'    - Keep taxa with ≥2 groups present in `group` and >4 samples total; drop
#'      `"unknown"`/`"uncultured"`.
#' 3. If `direction = "function"`:
#'    - Split by `function`; within each, map ASV to `taxon_level`, remove
#'      unclassified, and aggregate `norm_taxon_function_contrib` by level × sample.
#'    - Keep functions where `group` has ≥2 levels and >4 samples.
#' 4. If `pathway = "KO"`, strip `"ko:"` and convert KO abundances to KEGG pathway
#'    abundance (`ko2kegg_abundance`). For `"EC"`/`"MetaCyc"`, skip conversion.
#' 5. Build sample metadata subset from `sample_data(physeq)` for the samples present,
#'    then run `pathway_daa()` with the chosen `daa_method` and `p.adjust.method`.
#' 6. Annotate significant results via `MicroFun::pathway_annotation()` (for `"KO"`,
#'    pass `ko_to_kegg = TRUE`), or `result_annotation()` when `direction = "function"`.
#'
#' **Filtering rules:** taxa/functions are skipped if only one group is represented
#' in `group`, if total usable samples ≤ 4, or if the taxon label is
#' `"unknown"`/`"uncultured"`.
#'
#' @return
#' A named `list`. For `direction = "taxa"`, names are the merged taxon labels at
#' `taxon_level`; for `direction = "function"`, names are the function IDs being
#' tested. Each element is a `data.frame` of annotated significant features as
#' returned by `MicroFun::pathway_annotation()` or `result_annotation()`.
#' If no significant results are found, the function throws an error with a
#' suggestion to relax thresholds or check inputs.
#'
#' @section Expected input columns:
#' - `sample` (character): sample ID matching `sample_names(physeq)`.
#' - `taxon` (character): ASV/feature ID matching `taxa_names(physeq)`.
#' - `function` (character): KO/EC/MetaCyc identifier (e.g., `"ko:K00001"`, `"EC:1.1.1.1"`).
#' - `taxon_function_abun` (numeric): stratified abundance (for `direction = "taxa"`).
#' - `norm_taxon_function_contrib` (numeric): normalized contributions (for
#'   `direction = "function"`).
#'
#'
#' @import data.table reshape2 phyloseq
#' @export
#'
#' @examples
#' \dontrun{
#' # By taxon, KO → KEGG pathway analysis
#' res_ko <- pathway_sdaa(
#'   stratified_path = "picrust2_out/pathways_out/path_abun_strat.tsv.gz",
#'   physeq          = ps,
#'   taxon_level     = "Genus",
#'   group           = "condition",
#'   daa_method      = "Maaslin2",
#'   p.adjust.method = "BH",
#'   threshold       = 0.05,
#'   direction       = "taxa",
#'   pathway         = "KO"
#' )
#'
#' # By taxon, EC without conversion
#' res_ec <- pathway_sdaa(
#'   stratified_path = "picrust2_out/ec_out/ec_strat.tsv.gz",
#'   physeq          = ps,
#'   taxon_level     = "Genus",
#'   group           = "condition",
#'   direction       = "taxa",
#'   pathway         = "EC"
#' )
#'
#' # By function: aggregate ASVs to Genus within each function and test
#' res_func <- pathway_sdaa(
#'   stratified_path = "picrust2_out/pathways_out/path_abun_strat.tsv.gz",
#'   physeq          = ps,
#'   taxon_level     = "Genus",
#'   group           = "condition",
#'   direction       = "function",
#'   pathway         = "KO"
#' )
#' }
pathway_sdaa <- function(
    stratified_path,
    physeq,
    taxon_level = "Genus",
    group,
    daa_method = "Maaslin2",
    p.adjust.method = "none",
    threshold=0.05,
    direction = "taxa",
    pathway="KO"
) {
  # Load stratified Picrust2 output
  df <- data.table::fread(stratified_path, sep = "\t")

  sample_names <- phyloseq::sample_names(physeq)
  if (length(intersect(as.character(df$sample), phyloseq::sample_names(physeq))) == 0L)
  {
    stop("No overlapping samples between PICRUST2 output and phyloseq object.")
  }

  if (!("taxon" %in% names(df)) || !any(as.character(df$taxon) %in% phyloseq::taxa_names(physeq)))
  {
    stop("No overlapping taxa between PICRUSt2 'taxon' and phyloseq::taxa_names(physeq).")
  }
  mysam <- data.frame(physeq@sam_data)

  if (direction == "taxa") {

    # Split by taxon
    taxon_list <- split(df, df$taxon)

    # Convert each taxon into a wide-format table (KO x samples)
    result_list <- lapply(taxon_list, function(dt) {
      data.table::dcast(dt, `function` ~ sample, value.var = "norm_taxon_function_contrib", fill = 0)
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




