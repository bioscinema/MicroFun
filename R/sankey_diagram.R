#' Sankey diagram of genus -> pathway-class contributions (PICRUSt2 stratified)
#'
#' @description
#' Builds an interactive Sankey diagram (via **networkD3**) showing flows from
#' genera (left) to pathway superclasses (right), using per-genus-per-pathway
#' differential results (e.g., PICRUSt2 stratified output summarized by MaAsLin2).
#' Links are colored by effect direction (Up/Down). Node labels for left-side
#' genera can be rendered **to the left of the rectangles** for readability.
#'
#' @param pathway_sdaa_result A named \code{list} where each name is a genus and each
#'   element is a \code{data.frame} with at least the columns:
#'   \itemize{
#'     \item \code{pathway_class} Character; pathway superclass label (e.g., KEGG level-1).
#'     \item \code{logFC} Numeric; effect sign/size used to color links (>\,0 = Up).
#'     \item \code{p_adjust} Numeric; adjusted p-value (used only for filtering in user code).
#'   }
#'   Rows with \code{NA} in \code{pathway_class} are dropped.
#' @param taxa_list Optional character vector of genus names to include. Matching is
#'   case-insensitive. If supplied and none match, the function errors.
#' @param node_color Named character vector of length 2 giving link/node colors for
#'   \code{c(Up, Down)} directions. Default:
#'   \code{c(Up = "#ffa551", Down = "#70afdf")}.
#' @param node_neutral Character scalar used for neutral node color ("All"). Default \code{"#999999"}.
#'   (Links still use Up/Down colors.)
#' @param label_left_sources Logical; if \code{TRUE} (default), place genus labels on the
#'   **left** of left-column nodes and pathway labels on the **right** of right-column nodes.
#'   Implemented with a small \code{htmlwidgets::onRender()} D3 hook.
#' @param margin A named list with \code{top}, \code{right}, \code{bottom}, \code{left}
#'   margins (pixels) passed to \code{networkD3::sankeyNetwork()}. Increase \code{left}
#'   when \code{label_left_sources = TRUE} to avoid clipping.
#'
#' @details
#' The function constructs a long table of edges (genus -> pathway class), assigns a link
#' group \code{"Up"} or \code{"Down"} by the sign of \code{logFC}, builds unique node names,
#' and then calls \code{networkD3::sankeyNetwork()}. A post-render D3 callback adjusts text
#' anchors so that left-side labels appear to the left of their rectangles while right-side
#' labels remain on the right. Nodes are assigned a neutral group \code{"All"} by default,
#' so only links are colored by direction; you can customize node coloring later if desired.
#'
#' @return An \emph{htmlwidget} Sankey graph (class \code{htmlwidget}) that can be viewed
#'   in RStudio Viewer or a web browser, and embedded in R Markdown/Quarto.
#'
#'
#' @seealso \code{\link[networkD3]{sankeyNetwork}}, \code{\link[htmlwidgets]{onRender}}
#'
#' @importFrom networkD3 sankeyNetwork
#' @importFrom htmlwidgets JS onRender
#' @importFrom jsonlite toJSON
#' @export
tax2fun_sankey <- function(pathway_sdaa_result,
                           taxa_list = NULL,
                           node_color   = c(Up = "#ffa551", Down = "#70afdf"),
                           node_neutral = "#999999",
                           label_left_sources = TRUE,                 # <-- NEW
                           margin = list(top=10, right=140, bottom=10, left=160)) {

  # build long table (unchanged)
  sankey_df <- do.call(rbind, lapply(names(pathway_sdaa_result), function(genus) {
    df <- pathway_sdaa_result[[genus]]
    df <- df[!is.na(df$pathway_class), ]
    df$pathway_class <- sub(";.*", "", df$pathway_class)
    if (nrow(df) == 0) return(NULL)
    data.frame(
      source = rep(genus, nrow(df)),
      target = df$pathway_class,
      logFC  = df$logFC,
      qval   = df$p_adjust,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(sankey_df) || nrow(sankey_df) == 0) stop("No genus-function data.")

  if (!is.null(taxa_list)) {
    taxa_list <- unique(as.character(taxa_list))
    keep_mask <- tolower(sankey_df$source) %in% tolower(taxa_list)
    if (!any(keep_mask)) stop("None of the taxa in 'taxa_list' are present in the data.")
    sankey_df <- sankey_df[keep_mask, , drop = FALSE]
  }

  # nodes (+ side flag)
  nodes <- data.frame(name = unique(c(sankey_df$source, sankey_df$target)),
                      stringsAsFactors = FALSE)
  src_names     <- unique(sankey_df$source)
  nodes$side    <- ifelse(nodes$name %in% src_names, "left", "right")  # <-- NEW
  nodes$group   <- "All"

  # links
  sankey_df$group     <- ifelse(sankey_df$logFC > 0, "Up", "Down")
  sankey_df$source_id <- match(sankey_df$source, nodes$name) - 1
  sankey_df$target_id <- match(sankey_df$target, nodes$name) - 1
  sankey_df$value     <- 1

  # colors
  domain <- c("All","Up","Down")
  range  <- c(node_neutral, unname(node_color[c("Up","Down")]))
  color_scale <- htmlwidgets::JS(
    sprintf("d3.scaleOrdinal().domain(%s).range(%s)",
            jsonlite::toJSON(domain, auto_unbox=TRUE),
            jsonlite::toJSON(range,  auto_unbox=TRUE))
  )

  # base widget
  p <- networkD3::sankeyNetwork(
    Links  = sankey_df,
    Nodes  = nodes,
    Source = "source_id",
    Target = "target_id",
    Value  = "value",
    NodeID = "name",
    LinkGroup = "group",
    NodeGroup = "group",
    colourScale = color_scale,
    fontSize = 18,
    nodeWidth = 28,
    sinksRight = TRUE,
    margin = margin
  )

  # move labels: left nodes to left; right nodes to right
  if (label_left_sources) {
    p <- htmlwidgets::onRender(p, "
      function(el, x){
        var svg = d3.select(el).select('svg');
        svg.selectAll('.node text')
          .attr('text-anchor', function(d){
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? 'end' : 'start';
          })
          .attr('x', function(d){
            var w = (d.x1 !== undefined && d.x0 !== undefined) ? (d.x1 - d.x0)
                    : (d.dx !== undefined ? d.dx : 24);
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? -6 : (w + 6);
          });
      }
    ")
  }

  p
}


