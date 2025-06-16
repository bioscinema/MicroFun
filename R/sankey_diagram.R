#' Visualize Functional Contributions with a Sankey Diagram
#'
#' This function generates an interactive Sankey diagram to visualize the relationship
#' between genera and functional pathway classes, using the results from \code{pathway_sdaa()}.
#'
#' @param pathway_sdaa_result A named list of differential abundance results, typically the output from \code{pathway_sdaa()}.
#' @param node_color Named vector of colors to represent "Up" and "Down" regulated links. Default: \code{c(Up = "#79e8d0", Down = "#e8d079")}.
#'
#' @details
#' Each genus is treated as a source node and each pathway class as a target node.
#' The direction of regulation is determined by the sign of \code{logFC} in the input result.
#' Only the primary (first-level) pathway class is retained if a hierarchy exists.
#'
#' Internally, the function uses the \code{networkD3::sankeyNetwork()} function to render
#' an interactive diagram in the RStudio Viewer or web browser.
#'
#' @return A Sankey diagram object rendered by \code{networkD3::sankeyNetwork()}.
#'
#' @importFrom dplyr mutate
#' @importFrom networkD3 sankeyNetwork JS
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pathway_sdaa(...)
#' sankey_diagram(result)
#' }
sankey_diagram <- function(pathway_sdaa_result,
                           node_color = c(Up = "#79e8d0", Down = "#e8d079")) {
  # library(dplyr)
  # library(networkD3)

  # Flatten the result list into a single data frame
  sankey_df <- do.call(rbind, lapply(names(pathway_sdaa_result), function(genus) {
    df <- pathway_sdaa_result[[genus]]
    # 1) drop rows where pathway_class is NA
    df <- df[!is.na(df$pathway_class), ]
    # 2) keep only what's before the first “;”
    df$pathway_class <- sub(";.*", "", df$pathway_class)

    # if nothing remains, skip
    if (nrow(df) == 0) return(NULL)
    data.frame(
      source = rep(genus, nrow(df)),
      target = df$pathway_class,
      logFC = df$logFC,
      qval = df$p_adjust,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(sankey_df) || nrow(sankey_df) == 0) {
    stop("No genus-function associations found in input.")
  }

  # Create nodes and index
  nodes <- data.frame(name = unique(c(sankey_df$source, sankey_df$target)), stringsAsFactors = FALSE)
  nodes$group <- "All"  # Set all nodes to same group

  sankey_df <- sankey_df %>%
    mutate(source_id = match(source, nodes$name) - 1,
           target_id = match(target, nodes$name) - 1,
           value = 1,
           group = ifelse(logFC > 0, "Up", "Down"))

  # Define color scale JS
  color_scale <- JS("d3.scaleOrdinal()
    .domain(['All', 'Up', 'Down'])
    .range(['#999999', '#D73027', '#4575B4'])")

  # Plot Sankey diagram
  sankeyNetwork(Links = sankey_df,
                Nodes = nodes,
                Source = "source_id",
                Target = "target_id",
                Value = "value",
                NodeID = "name",
                LinkGroup = "group",
                NodeGroup = "group",
                colourScale = color_scale,
                fontSize = 20,
                nodeWidth = 50,
                sinksRight = FALSE)
}


