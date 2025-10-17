#' Dendrogram of top stratified functions per taxon (circular)
#'
#' @description
#' Build a **circular dendrogram** that connects each selected taxon (e.g., Genus)
#' to its **top-|logFC| functions** from a stratified pathway analysis result
#' (the list returned by `pathway_sdaa()` or a similar routine). Leaf nodes are
#' function labels colored by taxon, with point size proportional to
#' \eqn{|logFC|}. Ranking and selection use **absolute logFC only** (no p-value
#' filtering).
#'
#' @param pathway_sdaa_result A **named list** where each element corresponds to a
#'   taxon (e.g., Genus) and contains a `data.frame` of differential results with
#'   at least a `logFC` column and one label column among
#'   `pathway_name`, `feature`, or `function` (the function tries them in that order).
#' @param taxa_list Optional character vector of taxon names to include (must match
#'   names of `pathway_sdaa_result`). If `NULL`, all taxa in the list are used.
#' @param top_n_per_genus Integer. Number of functions to keep **per taxon** after
#'   ranking by `abs(logFC)`. Default `5`.
#' @param wrap_width Integer. Maximum character width for wrapping function labels.
#'   Default `36`.
#' @param title Character plot title. Default `"Function dendrogram"`.
#' @param point_size_range Numeric length-2 vector specifying the range of point
#'   sizes mapped to \eqn{|logFC|}. Default `c(1.5, 5)`.
#' @param label_size Numeric. Text size for leaf labels. Default `2`.
#'
#' @details
#' The function constructs a simple hierarchy: a root node → each selected taxon →
#' that taxon’s top functions. It then draws a circular dendrogram using **ggraph**
#' with edges colored by the **parent taxon** and leaf points sized by
#' \eqn{|logFC|}. Labels are wrapped to `wrap_width`.
#' **Note:** selection is based solely on `abs(logFC)`; if you want to enforce
#' significance, pre-filter `pathway_sdaa_result` before calling this function.
#'
#' @return
#' A named `list` with components:
#' \item{plot}{A `ggraph`/`ggplot` object (the circular dendrogram).}
#' \item{used_taxa}{Character vector of taxa actually plotted.}
#' \item{missing_taxa}{Character vector of requested taxa not found in `pathway_sdaa_result`.}
#'
#' @section Expected columns per taxon table:
#' \describe{
#'   \item{logFC}{Numeric signed log fold-change used for ranking and point size.}
#'   \item{pathway_name / feature / function}{One of these will be used as the label (in this order).}
#' }
#'
#' @seealso
#' `plot_pathway_logfc()` for barplots of signed effects; `ggraph::create_layout()`
#' for dendrogram layouts.
#'
#' @examples
#' \donttest{
#' # Suppose `res_list` is a named list of per-genus differential tables
#' # (each with columns like: pathway_name, feature, logFC).
#' # p_out <- plot_function_dendrogram(res_list, top_n_per_genus = 8)
#' # p_out$plot
#' }
#'
#' @export
plot_function_dendrogram <- function(pathway_sdaa_result,
                                           taxa_list       = NULL,
                                           top_n_per_genus = 5,
                                           wrap_width      = 36,
                                           title           = "Function dendrogram",
                                           point_size_range= c(1.5, 5),
                                           label_size      = 2) {
  label_col       = c("pathway_name","feature","function")
  palette         = "Paired"
  stopifnot(length(pathway_sdaa_result) > 0)
  label_col <- label_col[1]

  # ---- filter by taxa_list if provided ----
  all_genera <- names(pathway_sdaa_result)
  if (!is.null(taxa_list)) {
    taxa_list <- unique(as.character(taxa_list))
    used_genera    <- intersect(taxa_list, all_genera)
    missing_genera <- setdiff(taxa_list, all_genera)
    if (!length(used_genera)) stop("None of the taxa in `taxa_list` are found in `pathway_sdaa_result`.")
    pathway_sdaa_result <- pathway_sdaa_result[used_genera]
  } else {
    used_genera    <- all_genera
    missing_genera <- character(0)
  }

  # ---- collect rows; rank ONLY by |logFC|; NO padj anywhere ----
  df_all <- dplyr::bind_rows(lapply(names(pathway_sdaa_result), function(gn) {
    df <- pathway_sdaa_result[[gn]]
    if (is.null(df) || !nrow(df)) return(NULL)
    lbl <- df[[label_col]]; if (is.null(lbl)) lbl <- df[["feature"]]
    dplyr::tibble(
      genus     = gn,
      label_raw = lbl,
      logFC     = df$logFC
    )
  })) %>%
    dplyr::filter(!is.na(logFC)) %>%
    dplyr::group_by(genus) %>%
    dplyr::arrange(dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_per_genus) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label    = stringr::str_wrap(dplyr::if_else(is.na(label_raw) | label_raw == "", "(unknown)", label_raw),
                                   width = wrap_width),
      size_val = abs(logFC)
    )

  if (!nrow(df_all)) return(list(plot=NULL, used_taxa=used_genera, missing_taxa=missing_genera))

  # ---- edges & vertices ----
  root <- "origin"
  edges_g <- dplyr::distinct(df_all, from = root, to = genus)
  edges_f <- df_all %>% dplyr::transmute(from = genus, to = paste(genus, label, sep = " :: "))
  edges   <- dplyr::bind_rows(edges_g, edges_f)

  genus_levels <- unique(df_all$genus)
  leaf_names   <- unique(edges_f$to)

  vertices <- dplyr::tibble(name = unique(c(edges$from, edges$to))) %>%
    dplyr::mutate(
      is_genus = name %in% genus_levels,
      is_leaf  = name %in% leaf_names,
      genus    = dplyr::case_when(
        name %in% genus_levels ~ name,
        name %in% leaf_names   ~ sub(" ::.*$", "", name),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::left_join(
      df_all %>% dplyr::mutate(node = paste(genus, label, sep = " :: ")) %>%
        dplyr::select(node, logFC, size_val),
      by = c("name" = "node")
    )

  # ---- colors ----
  n_gen <- length(genus_levels)
  if (n_gen <= 12 && palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    gen_cols <- RColorBrewer::brewer.pal(max(3, min(12, n_gen)), palette)[seq_len(n_gen)]
  } else {
    gen_cols <- scales::hue_pal()(n_gen)
  }
  names(gen_cols) <- genus_levels

  # ---- graph & layout ----
  g <- igraph::graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
  em <- igraph::ends(g, igraph::E(g), names = TRUE)
  igraph::E(g)$branch_col <- ifelse(em[,1] == root, "white", gen_cols[ em[,1] ])

  V(g)$is_leaf  <- vertices$is_leaf
  V(g)$genus    <- vertices$genus
  V(g)$logFC    <- vertices$logFC
  V(g)$size_val <- vertices$size_val

  layout_df <- ggraph::create_layout(g, layout = "dendrogram", circular = TRUE) %>%
    dplyr::mutate(
      angle = 180 / pi * atan2(y, x),
      angle = ifelse(angle < 0, angle + 360, angle),
      hjust = ifelse(angle > 90 & angle < 270, 1, 0),
      angle = ifelse(angle > 90 & angle < 270, angle + 180, angle)
    )

  # ---- plot ----
  p <- ggraph::ggraph(layout_df) +
    ggraph::geom_edge_diagonal(aes(colour = I(igraph::E(g)$branch_col)),
                               width = 0.8, alpha = 0.7) +
    ggraph::geom_node_point(aes(filter = !is_leaf & !is.na(genus),
                                x = x*1.00, y = y*1.00,
                                colour = genus),
                            size = 0.4, alpha = 0.9) +
    ggraph::geom_node_point(aes(filter = is_leaf,
                                x = x*1.06, y = y*1.06,
                                colour = genus,
                                size = size_val),
                            alpha = 0.7, show.legend = FALSE) +
    ggraph::geom_node_text(aes(filter = is_leaf,
                               x = x*1.12, y = y*1.12,
                               label = sub("^.* :: ", "", name),
                               angle = angle, hjust = hjust,
                               colour = genus),
                           size = label_size, alpha = 0.95, show.legend = FALSE) +
    scale_colour_manual(values = gen_cols, name = "Genus",
                        guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    scale_size_continuous(range = point_size_range, guide = "none") +  # <-- no size legend
    coord_equal() +
    theme_void() +
    theme(legend.position = "right",
          plot.title.position = "plot",
          plot.margin = unit(c(5,5,5,5), "mm")) +
    ggtitle(title) +
    expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))

  list(plot = p, used_taxa = used_genera, missing_taxa = missing_genera)
}


