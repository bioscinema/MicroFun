#' Plot signed logFC pathways for each taxon (list of plots)
#'
#' @description
#' Given the list returned by `pathway_sdaa()` (one data frame per taxon/function
#' group), make a bar plot of signed log fold-changes (logFC) for each element,
#' highlighting the direction of change and labeling bars with wrapped pathway
#' names. Returns a named list of `ggplot` objects.
#'
#' @param pathway_sdaa_result A **named list** where each element is a
#'   `data.frame` of differential results for one taxon (or function), containing
#'   at least the columns `feature`, `logFC`, `p_adjust`, `group1`, `group2`,
#'   and optionally `pathway_name`.
#' @param top_n Integer. Maximum number of bars (features) to display per plot.
#'   The function keeps the strongest absolute effects, but always keeps
#'   features with `p_adjust <= padj_cutoff`. Default `15`.
#' @param padj_cutoff Numeric. Adjusted p-value cutoff for retaining significant
#'   features even if they are not among the top `top_n` by absolute effect.
#'   Default `0.05`.
#' @param wrap_width Integer. Character width used to line-wrap labels (via
#'   `stringr::str_wrap`). Default `45`.
#'
#' @return A named `list` of `ggplot` objects, one per entry of
#'   `pathway_sdaa_result`. Entries with missing or empty data are omitted.
#'
#' @details
#' Internally calls [plot_pathway_logfc_one()] on each element using
#' `label_col = "pathway_name"`; when that column is missing or empty the
#' function falls back to `feature`.
#'
#' @seealso [plot_pathway_logfc_one()]
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip geom_hline scale_fill_manual labs theme_minimal theme element_text margin
#' @importFrom dplyr mutate arrange filter slice_head
#' @importFrom forcats fct_reorder
#' @importFrom stringr str_wrap
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
plot_pathway_logfc <- function(pathway_sdaa_result,
                               top_n = 15,
                               padj_cutoff = 0.05,
                               wrap_width = 45) {
  taxa = names(pathway_sdaa_result)
  label_col = "pathway_name"
  out <- lapply(taxa, function(tn) {
    if (!tn %in% names(pathway_sdaa_result)) return(NULL)
    plot_pathway_logfc_one(
      df = pathway_sdaa_result[[tn]],
      title = tn,
      top_n = top_n,
      padj_cutoff = padj_cutoff,
      label_col = label_col,
      wrap_width = wrap_width
    )
  })
  names(out) <- taxa
  Filter(Negate(is.null), out)
}


plot_pathway_logfc_one <- function(df,
                                   title = NULL,
                                   padj_cutoff = 0.05,
                                   top_n = top_n,
                                   label_col = c("pathway_name", "feature"),
                                   wrap_width = 45) {
  label_col <- label_col[1]

  df2 <- df %>%
    mutate(
      p_adjust = ifelse(is.na(p_adjust), 1, p_adjust),
      label = .data[[label_col]],
      label = ifelse(is.na(label) | label == "", .data[["feature"]], label),
      label = stringr::str_wrap(label, width = wrap_width)
    ) %>%
    arrange(desc(abs(logFC))) %>%
    # keep the strongest effects, but always keep significant ones
    filter(row_number() <= top_n | p_adjust <= padj_cutoff) %>%
    slice_head(n = top_n) %>%
    filter(!is.na(logFC)) %>%
    mutate(label = fct_reorder(label, logFC))

  g1 <- ggplot(df2, aes(x = label, y = logFC, fill = logFC > 0)) +
    geom_col(width = 0.8) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("FALSE" = "#70afdf", "TRUE" = "#ffa551"), guide = "none") +
    labs(
      x = NULL, y = "logFC",
      title = title %||% "",
      subtitle = paste0(unique(df$group1), " vs ", unique(df$group2))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.margin = margin(10, 18, 10, 36),   # add left margin to avoid clipping
      axis.text.y = element_text(hjust = 1)
    )
  g1
}
