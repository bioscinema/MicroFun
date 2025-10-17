#' Circular dendrogram of top-N functions per taxon
#'
#' Create a circular dendrogram where each genus (or other taxon label)
#' branches to its top-N associated functions (ranked by |logFC|). Node
#' point size is proportional to |logFC|, and labels can be drawn from
#' different columns (pathway name, feature ID, or function description).
#'
#' @param pathway_sdaa_result Named list of results from \code{pathway_sdaa()}.
#'   Each element should be a data.frame with at least a \code{logFC} column
#'   and one of \code{pathway_name}, \code{feature}, or \code{function}.
#'   The list names are used as genera/taxa labels.
#' @param taxa_list Optional character vector of genera to include. If \code{NULL}
#'   (default), all taxa in \code{pathway_sdaa_result} are included.
#'   If provided, missing taxa are reported in the return object.
#' @param top_n_per_genus Integer. Number of functions per genus to display
#'   (ranked by absolute \code{logFC}). Default: 5.
#' @param wrap_width Integer. Width for wrapping long function labels.
#'   Passed to \code{stringr::str_wrap()}. Default: 36.
#' @param title Character string, plot title. Default: \code{"Function dendrogram"}.
#' @param point_size_range Numeric length-2 vector. Range of point sizes for leaf
#'   nodes, scaled to |logFC| values. Default: \code{c(1.5, 5)}.
#' @param label_size Numeric. Text size for function labels. Default: 2.
#' @param label_priority Character. Which column to use for labeling functions.
#'   Must be one of \code{"pathway_name"}, \code{"feature"}, or \code{"function"}.
#'   If the chosen column is absent, the function falls back to \code{feature}
#'   or \code{function}. Default: \code{"pathway_name"}.
#' @param palette Character. Name of RColorBrewer palette to use for genus colors
#'   (if \eqn{\leq 12} genera). Defaults to \code{"Paired"}. If there are more
#'   than 12 genera, falls back to a hue-based palette from \pkg{scales}.
#'
#' @details
#' For each genus in the input, the function selects the top-N functions
#' ranked by absolute log fold change (\code{|logFC|}). It builds a dendrogram
#' with a single root node branching into genera nodes, which then branch into
#' leaf nodes (functions). Each leaf node is plotted with a size proportional
#' to \code{|logFC|}, and colored according to its genus.
#'
#' @return A list with three elements:
#'   \item{plot}{A \code{ggplot} dendrogram object (from \pkg{ggraph}).}
#'   \item{used_taxa}{Character vector of genera included in the plot.}
#'   \item{missing_taxa}{Character vector of requested taxa (via \code{taxa_list})
#'   that were not found in \code{pathway_sdaa_result}.}
#'
#' @examples
#' \dontrun{
#' # Fake input
#' res <- list(
#'   GenusA = data.frame(pathway_name = c("Glycolysis","TCA cycle"),
#'                       logFC = c(1.2,-0.8)),
#'   GenusB = data.frame(feature = c("K00001","K00002"),
#'                       logFC = c(2.1, 0.5))
#' )
#'
#' dend <- tax_fun_dendrogram(res, top_n_per_genus = 2,
#'                            label_priority = "pathway_name")
#' dend$plot
#' }
#'
#' @seealso \code{\link{pathway_sdaa}}, \pkg{ggraph}, \pkg{igraph}
#'
#' @export
tax_fun_dendrogram <- function(
    pathway_sdaa_result,
    taxa_list        = NULL,
    top_n_per_genus  = 5,
    wrap_width       = 36,
    title            = "Function dendrogram",
    point_size_range = c(1.5, 5),
    label_size       = 2,
    label_priority   = "pathway_name",
    palette          = "Paired"
){
  # ---- basic checks ------------------------------------------------------------
  stopifnot(is.list(pathway_sdaa_result), length(pathway_sdaa_result) > 0)
  stopifnot(is.numeric(top_n_per_genus), length(top_n_per_genus) == 1, top_n_per_genus > 0)
  stopifnot(is.numeric(point_size_range), length(point_size_range) == 2)
  label_priority <- match.arg(label_priority)

  all_genera <- names(pathway_sdaa_result)
  if (is.null(all_genera) || anyNA(all_genera) || any(all_genera == "")) {
    stop("`pathway_sdaa_result` must be a named list (names = genera).")
  }

  # Subset by taxa_list (preserve user-given order)
  if (!is.null(taxa_list)) {
    taxa_list <- unique(as.character(taxa_list))
    used_genera    <- intersect(taxa_list, all_genera)
    missing_genera <- setdiff(taxa_list, all_genera)
    if (!length(used_genera)) stop("None of `taxa_list` are in `pathway_sdaa_result`.")
    pathway_sdaa_result <- pathway_sdaa_result[used_genera]
  } else {
    used_genera    <- all_genera
    missing_genera <- character(0)
  }

  # Helper: resolve the label column with fallbacks
  pick_label <- function(df){
    if (!is.data.frame(df) || !nrow(df)) return(rep(NA_character_, 0))
    cols <- c(label_priority, "feature", "function")  # ensure sensible fallback
    cols <- intersect(cols, colnames(df))
    if (!length(cols)) return(rep(NA_character_, nrow(df)))
    df[[cols[1]]]
  }

  # ---- collect & rank rows (top-N by |logFC| per genus) ------------------------
  df_list <- lapply(names(pathway_sdaa_result), function(gn){
    df <- pathway_sdaa_result[[gn]]
    if (is.null(df) || !is.data.frame(df) || !"logFC" %in% names(df) || !nrow(df)) return(NULL)
    tibble::tibble(
      genus     = gn,
      label_raw = pick_label(df),
      logFC     = df$logFC
    )
  })
  df_all <- dplyr::bind_rows(df_list)

  if (is.null(df_all) || !nrow(df_all)) {
    return(list(plot=NULL, used_taxa=used_genera, missing_taxa=missing_genera))
  }

  df_all <- df_all |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::mutate(
      label_raw = dplyr::if_else(is.na(label_raw) | label_raw == "", "(unknown)", label_raw),
      genus     = factor(genus, levels = used_genera)   # keep input order on the ring
    ) |>
    dplyr::group_by(genus) |>
    dplyr::arrange(dplyr::desc(abs(logFC)), .by_group = TRUE) |>
    dplyr::slice_head(n = top_n_per_genus) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      label    = stringr::str_wrap(label_raw, width = wrap_width),
      size_val = abs(logFC)
    )

  if (!nrow(df_all)) {
    return(list(plot=NULL, used_taxa=used_genera, missing_taxa=missing_genera))
  }

  # ---- build tree edges & vertices --------------------------------------------
  root    <- "origin"
  edges_g <- dplyr::distinct(df_all, from = root, to = genus)
  edges_f <- df_all |>
    dplyr::transmute(from = as.character(genus), to = paste(genus, label, sep = " :: "))

  edges <- dplyr::bind_rows(edges_g, edges_f)

  genus_levels <- levels(df_all$genus)
  leaf_names   <- unique(edges_f$to)

  vertices <- tibble::tibble(name = unique(c(edges$from, edges$to))) |>
    dplyr::mutate(
      is_genus = name %in% genus_levels,
      is_leaf  = name %in% leaf_names,
      genus    = dplyr::case_when(
        name %in% genus_levels ~ name,
        name %in% leaf_names   ~ sub(" ::.*$", "", name),
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::left_join(
      df_all |>
        dplyr::mutate(node = paste(genus, label, sep = " :: ")) |>
        dplyr::select(node, logFC, size_val),
      by = c("name" = "node")
    )

  # ---- colors ------------------------------------------------------------------
  n_gen <- length(genus_levels)
  if (n_gen <= 12 && palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    gen_cols <- RColorBrewer::brewer.pal(max(3, min(12, n_gen)), palette)[seq_len(n_gen)]
  } else {
    gen_cols <- scales::hue_pal()(n_gen)
  }
  names(gen_cols) <- genus_levels

  # ---- graph + attributes ------------------------------------------------------
  g <- igraph::graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

  # branch color: edges from root are white, otherwise color by parent (genus node)
  em <- igraph::ends(g, igraph::E(g), names = TRUE)
  parent_name <- em[,1]
  branch_col  <- ifelse(parent_name == root, "white", gen_cols[parent_name])
  igraph::E(g)$branch_col <- branch_col

  # copy vertex attrs
  igraph::V(g)$is_leaf  <- vertices$is_leaf
  igraph::V(g)$genus    <- vertices$genus
  igraph::V(g)$logFC    <- vertices$logFC
  igraph::V(g)$size_val <- vertices$size_val

  # ---- layout ------------------------------------------------------------------
  layout_df <- ggraph::create_layout(g, layout = "dendrogram", circular = TRUE) |>
    dplyr::mutate(
      angle = (180 / pi) * atan2(y, x),
      angle = ifelse(angle < 0, angle + 360, angle),
      hjust = ifelse(angle > 90 & angle < 270, 1, 0),
      angle = ifelse(angle > 90 & angle < 270, angle + 180, angle)
    )

  # ---- plot --------------------------------------------------------------------
  p <- ggraph::ggraph(layout_df) +
    ggraph::geom_edge_diagonal(aes(colour = I(igraph::E(g)$branch_col)),
                               width = 0.8, alpha = 0.7) +
    ggraph::geom_node_point(
      aes(filter = !is_leaf & !is.na(genus), x = x*1.00, y = y*1.00, colour = genus),
      size = 0.6, alpha = 0.9
    ) +
    ggraph::geom_node_point(
      aes(filter = is_leaf, x = x*1.06, y = y*1.06, colour = genus, size = size_val),
      alpha = 0.7, show.legend = FALSE
    ) +
    ggraph::geom_node_text(
      aes(filter = is_leaf, x = x*1.12, y = y*1.12,
          label = sub("^.* :: ", "", name), angle = angle, hjust = hjust, colour = genus),
      size = label_size, alpha = 0.95, show.legend = FALSE, lineheight = 0.9
    ) +
    ggplot2::scale_colour_manual(values = gen_cols, name = "Genus",
                                 guide = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))) +
    ggplot2::scale_size_continuous(range = point_size_range, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title.position = "plot",
      plot.margin = grid::unit(c(5,5,5,5), "mm")
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))

  list(plot = p, used_taxa = used_genera, missing_taxa = missing_genera)
}


