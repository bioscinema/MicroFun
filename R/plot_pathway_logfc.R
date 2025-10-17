#' Plot logFC barplots for pathway_sdaa results
#'
#' Generate per-taxon barplots showing log fold changes (logFC) of the top-N
#' pathways/functions returned by \code{pathway_sdaa()}. Each plot is a
#' horizontal bar chart, optionally filtering by adjusted p-values and wrapping
#' long labels.
#'
#' @param pathway_sdaa_result Named list of results from \code{pathway_sdaa()}.
#'   Each element should be a data.frame with at least a \code{logFC} column
#'   and a feature label column (e.g., \code{pathway_name}, \code{feature}).
#'   List names are used as taxa/taxon names.
#' @param top_n Integer. Maximum number of features (per taxon) to plot.
#'   Ranked by absolute logFC. Default: 15.
#' @param padj_cutoff Numeric. Adjusted p-value cutoff; features passing this
#'   threshold are always retained even if outside the top-N. Default: 0.05.
#' @param wrap_width Integer. Width used for wrapping long labels with
#'   \code{stringr::str_wrap()}. Default: 45.
#' @param label_col Character. Preferred label column to display on the y-axis.
#'   Options: \code{"pathway_name"} or \code{"feature"}. Default:
#'   \code{"pathway_name"}.
#' @param pos_color,neg_color Character. Fill colors for positive vs negative
#'   logFC bars. Defaults: orange (\code{"#ffa551"}) and blue
#'   (\code{"#70afdf"}).
#'
#'
#' @return
#' \code{tax_fun_logfc()} returns a named list of \code{ggplot} objects,
#' one per taxon. \code{plot_pathway_logfc_one()} returns a single \code{ggplot}.
#'
#' @examples
#' \dontrun{
#' # Example input (toy data)
#' res <- list(
#'   GenusA = data.frame(pathway_name = c("Glycolysis","TCA cycle"),
#'                       feature = c("K001","K002"),
#'                       logFC = c(1.2,-0.8),
#'                       p_adjust = c(0.01, 0.2),
#'                       group1 = "Case", group2 = "Control"),
#'   GenusB = data.frame(feature = c("K01001","K01002"),
#'                       logFC = c(-1.5, 0.7),
#'                       p_adjust = c(0.03, 0.8),
#'                       group1 = "Case", group2 = "Control")
#' )
#'
#' plots <- tax_fun_logfc(res, top_n = 5)
#' plots$GenusA   # ggplot object for GenusA
#' }
#'
#' @seealso \code{\link{pathway_sdaa}}, \code{ggplot2}
#'
#' @export
tax_fun_logfc <- function(
    pathway_sdaa_result,
    top_n       = 15,
    padj_cutoff = 0.05,
    wrap_width  = 45,
    label_col   = c("pathway_name", "feature"),   # preferred label column(s)
    pos_color   = "#ffa551",
    neg_color   = "#70afdf"
){
  stopifnot(is.list(pathway_sdaa_result), length(pathway_sdaa_result) > 0)
  taxa <- names(pathway_sdaa_result)
  if (is.null(taxa) || anyNA(taxa) || any(taxa == "")) {
    stop("`pathway_sdaa_result` must be a *named* list (names = taxa).")
  }

  label_col <- match.arg(label_col)

  out <- lapply(taxa, function(tn){
    df <- pathway_sdaa_result[[tn]]
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)

    plot_pathway_logfc_one(
      df          = df,
      title       = tn,
      top_n       = top_n,
      padj_cutoff = padj_cutoff,
      label_col   = label_col,
      wrap_width  = wrap_width,
      pos_color   = pos_color,
      neg_color   = neg_color
    )
  })
  names(out) <- taxa
  Filter(Negate(is.null), out)
}

# ---- single-panel helper -------------------------------------------------------
plot_pathway_logfc_one <- function(
    df,
    title        = NULL,
    padj_cutoff  = 0.05,
    top_n        = 15,
    label_col    = c("pathway_name", "feature"),
    wrap_width   = 45,
    pos_color    = "#ffa551",
    neg_color    = "#70afdf"
){
  label_col <- match.arg(label_col)

  # required columns
  req <- c("logFC")
  if (!all(req %in% names(df))) {
    stop("`df` must contain: ", paste(req, collapse = ", "))
  }

  # locate an adjusted-p column (best-effort)
  padj_candidates <- intersect(c("p_adjust","padj","FDR","qvalue","q_value","qval"), names(df))
  padj_col <- if (length(padj_candidates)) padj_candidates[1] else NA_character_

  # group labels (optional, for subtitle)
  g1 <- if ("group1" %in% names(df)) unique(df$group1) else NA
  g2 <- if ("group2" %in% names(df)) unique(df$group2) else NA
  subtitle_txt <- if (!anyNA(c(g1,g2))) paste0(g1, " vs ", g2) else NULL

  # pick label column, with fallback to 'feature'
  lbl <- df[[label_col]]
  if (is.null(lbl) || all(is.na(lbl) | lbl == "")) {
    lbl <- df[["feature"]]
  }
  if (is.null(lbl)) lbl <- rep("(unknown)", nrow(df))

  df2 <- dplyr::mutate(
    df,
    label_raw = lbl,
    p_adjust  = if (!is.na(padj_col)) .data[[padj_col]] else NA_real_,
    p_adjust  = dplyr::if_else(is.na(p_adjust), 1, p_adjust),
    label     = stringr::str_wrap(dplyr::if_else(is.na(label_raw) | label_raw == "", "(unknown)", label_raw),
                                  width = wrap_width)
  )

  # keep strongest effects, but always keep significant ones
  df2 <- df2 |>
    dplyr::arrange(dplyr::desc(abs(.data$logFC))) |>
    dplyr::filter(dplyr::row_number() <= top_n | .data$p_adjust <= padj_cutoff) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::filter(!is.na(.data$logFC)) |>
    dplyr::mutate(label = forcats::fct_reorder(.data$label, .data$logFC))

  if (!nrow(df2)) return(NULL)

  # sign fill as factor for stable legend mapping (even though legend is hidden)
  df2$sign <- factor(df2$logFC > 0, levels = c(FALSE, TRUE))

  ggplot2::ggplot(df2, ggplot2::aes(x = .data$label, y = .data$logFC, fill = .data$sign)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("FALSE" = neg_color, "TRUE" = pos_color), guide = "none") +
    ggplot2::labs(
      x = NULL, y = "logFC",
      title    = if (is.null(title)) "" else title,
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(10, 18, 10, 36),
      axis.text.y = ggplot2::element_text(hjust = 1)
    )
}

