#' Sankey diagram of functions (KO/EC/MetaCyc) to taxa
#'
#' Build an interactive Sankey diagram showing connections between functions
#' (KO, EC, or MetaCyc identifiers, optionally grouped to higher-level categories)
#' and taxa/features from differential abundance analysis results. Edges are
#' colored by the sign of log fold change (Up/Down), while nodes are colored
#' with neutral or user-specified palettes.
#'
#' @param pathway_sdaa_result Named list of results from \code{pathway_sdaa()}.
#'   Each element should be a data.frame with at least the columns
#'   \code{feature}, \code{logFC}, and an adjusted p-value column
#'   (\code{p_adjust}, \code{padj}, \code{FDR}, \code{qval}, etc.).
#'   List names must be function IDs (KO, EC, or MetaCyc identifiers).
#' @param pathway Character scalar. Functional namespace of the function IDs;
#'   one of \code{"KO"}, \code{"EC"}, or \code{"MetaCyc"}.
#' @param group_functions Logical. If \code{TRUE}, function IDs are mapped to
#'   higher-level group labels using reference mappings shipped in the package:
#'   KO to PathwayL1, EC/MetaCyc to description. Default: \code{TRUE}.
#' @param taxa_list Optional character vector. If supplied, only edges whose
#'   source (function group or raw function ID) matches one of these values
#'   will be retained in the diagram.
#' @param node_color Named character vector of length 2 giving fill colors for
#'   edges with positive logFC (\code{"Up"}) and negative logFC (\code{"Down"}).
#'   Default: orange for Up (\code{"#ffa551"}), blue for Down (\code{"#70afdf"}).
#' @param node_neutral Character. Color used for neutral nodes (the default
#'   group "All"). Default: \code{"#999999"}.
#' @param label_left_sources Logical. If \code{TRUE}, left-side node labels are
#'   right-aligned and wrapped to avoid overlap. Default: \code{TRUE}.
#' @param margin List of integers. Margins for the sankey plot in pixels, with
#'   elements \code{top}, \code{right}, \code{bottom}, \code{left}.
#'   Default: \code{list(top=10, right=200, bottom=10, left=180)}.
#' @param height,width Numeric. Output height and width in pixels for the
#'   interactive sankey widget.
#'
#' @examples
#' \dontrun{
#' # Example (toy data)
#' res <- list(
#'   K00001 = data.frame(feature = c("GenusA","GenusB"),
#'                       logFC   = c(1.2,-0.8),
#'                       p_adjust = c(0.01,0.2)),
#'   K00002 = data.frame(feature = c("GenusC"),
#'                       logFC   = c(-1.5),
#'                       p_adjust = c(0.04))
#' )
#'
#' p <- fun2tax_sankey(res, pathway = "KO", group_functions = FALSE)
#' p   # interactive sankey widget
#' }
#'
#' @seealso \code{\link{pathway_sdaa}}, \pkg{networkD3}, \pkg{htmlwidgets}
#'
#' @importFrom stats setNames
#' @importFrom jsonlite toJSON
#' @importFrom htmlwidgets JS onRender
#' @importFrom networkD3 sankeyNetwork
#' @export
fun2tax_sankey <- function(
    pathway_sdaa_result,                       # named list: names = function IDs; each df has feature, logFC, p_adjust
    pathway            = c("KO","EC","MetaCyc"),
    group_functions    = TRUE,                 # if TRUE, map function IDs to higher-level groups (KO: PathwayL1; EC/MetaCyc: description)
    taxa_list          = NULL,                 # optional filter: keep only edges whose SOURCE (group or function) matches
    node_color         = c(Up = "#ffa551", Down = "#70afdf"),
    node_neutral       = "#999999",
    label_left_sources = TRUE,
    margin             = list(top=10, right=200, bottom=10, left=180),
    height             = 1200,
    width              = 1000
){
  # ---- validations -------------------------------------------------------------
  if (!is.list(pathway_sdaa_result) || !length(pathway_sdaa_result)) {
    stop("`pathway_sdaa_result` must be a non-empty list.")
  }
  pathway <- match.arg(pathway)

  fun_ids_raw <- names(pathway_sdaa_result)
  if (is.null(fun_ids_raw) || anyNA(fun_ids_raw) || any(fun_ids_raw == "")) {
    stop("`pathway_sdaa_result` must be a *named* list (names = function IDs).")
  }

  if (!is.character(node_neutral) || length(node_neutral) != 1L) {
    stop("`node_neutral` must be a single color string.")
  }
  if (!is.character(node_color) || !all(c("Up","Down") %in% names(node_color))) {
    stop("`node_color` must be a named character vector with names 'Up' and 'Down'.")
  }

  # ---- helpers -----------------------------------------------------------------
  normalize_fun_ids <- function(ids, pathway = c("KO","EC","MetaCyc")) {
    pathway <- match.arg(pathway)
    ids <- trimws(as.character(ids))
    if (pathway == "KO") {
      ids <- toupper(gsub("^ko:", "", ids, ignore.case = TRUE))
    } else if (pathway == "EC") {
      ids <- gsub("^ec:", "", ids, ignore.case = TRUE)
      ids <- gsub("^EC[\\s_]*", "", ids, ignore.case = TRUE)
    } else { # MetaCyc
      ids <- gsub("^metacyc:", "", ids, ignore.case = TRUE)
    }
    ids
  }

  load_reference_map <- function(pathway) {
    file_name <- switch(pathway,
                        KO      = "KO_reference.RData",
                        EC      = "EC_reference.RData",
                        MetaCyc = "MetaCyc_reference.RData")
    ref_path <- system.file("extdata", file_name, package = "MicroFun")
    if (!nzchar(ref_path) || !file.exists(ref_path)) {
      stop("Reference file not found in MicroFun/extdata: ", file_name)
    }
    env <- new.env(parent = emptyenv())
    load(ref_path, envir = env)
    ref <- NULL
    for (nm in ls(env)) {
      cand <- get(nm, envir = env)
      if (is.data.frame(cand)) { ref <- cand; break }
    }
    if (is.null(ref)) stop("No data.frame found inside ", file_name)

    if (pathway == "KO") {
      if (!all(c("KO","PathwayL1") %in% names(ref))) stop("KO reference must have 'KO' and 'PathwayL1' columns.")
      stats::setNames(as.character(ref$PathwayL1), normalize_fun_ids(ref$KO, "KO"))
    } else if (pathway == "EC") {
      if (!all(c("EC","description") %in% names(ref))) stop("EC reference must have 'EC' and 'description' columns.")
      stats::setNames(as.character(ref$description), normalize_fun_ids(ref$EC, "EC"))
    } else {
      id_col <- if ("MetaCyc" %in% names(ref)) "MetaCyc" else if ("ID" %in% names(ref)) "ID" else NA_character_
      if (is.na(id_col) || !("description" %in% names(ref))) stop("MetaCyc reference must have ID and 'description' columns.")
      ids_norm <- normalize_fun_ids(ref[[id_col]], "MetaCyc")
      stats::setNames(as.character(ref$description), ids_norm)
    }
  }

  get_qcol <- function(df) {
    cand <- c("p_adjust","padj","qval","q_value","qvalue","p.adj","FDR","p_adj")
    hit <- intersect(cand, names(df))
    if (length(hit)) hit[1] else NULL
  }

  # ---- normalize function IDs on input list ------------------------------------
  fun_ids_norm <- normalize_fun_ids(fun_ids_raw, pathway)
  keep <- !is.na(fun_ids_norm) & nzchar(fun_ids_norm)
  if (!any(keep)) stop("After normalization, no valid function IDs remain.")
  pathway_sdaa_result <- pathway_sdaa_result[keep]
  fun_ids_norm <- fun_ids_norm[keep]
  names(pathway_sdaa_result) <- fun_ids_norm

  # ---- optional mapping from function IDs to higher-level groups ---------------
  id2group <- NULL
  if (isTRUE(group_functions)) {
    id2group <- load_reference_map(pathway)
  }

  # ---- assemble edges: source -> target (group/function -> taxon) --------------
  build_rows_from_df <- function(df) {
    # Return NULL if structure is unsuitable
    qcol <- get_qcol(df)
    if (is.null(qcol) || !("feature" %in% names(df)) || !("logFC" %in% names(df))) return(NULL)
    keep <- !is.na(df$feature) & nzchar(df$feature) & !is.na(df[[qcol]])
    if (!any(keep)) return(NULL)
    data.frame(
      feature  = as.character(df$feature[keep]),
      logFC    = as.numeric(df$logFC[keep]),
      p_adjust = as.numeric(df[[qcol]][keep]),
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(group_functions)) {
    # Map each function to a group label and aggregate within group
    grp_lab <- unname(id2group[fun_ids_norm])
    valid_fun <- !is.na(grp_lab) & nzchar(grp_lab)
    fun_ids_norm <- fun_ids_norm[valid_fun]
    grp_lab <- grp_lab[valid_fun]
    if (!length(fun_ids_norm)) stop("No function IDs matched the reference mapping after normalization.")

    groups <- split(fun_ids_norm, grp_lab)

    pieces <- lapply(names(groups), function(glab) {
      fids <- groups[[glab]]
      rows_list <- lapply(fids, function(fid) build_rows_from_df(pathway_sdaa_result[[fid]]))
      rows_list <- Filter(Negate(is.null), rows_list)
      if (!length(rows_list)) return(NULL)
      rows <- do.call(rbind, rows_list)
      if (!nrow(rows)) return(NULL)

      # Deduplicate per taxon within group:
      # keep the row with the LARGER p_adjust (matches original rule).
      ord <- order(rows$feature, rows$p_adjust, decreasing = TRUE)
      rows <- rows[ord, , drop = FALSE]
      rows <- rows[!duplicated(rows$feature), , drop = FALSE]

      data.frame(
        source = rep(glab, nrow(rows)),         # group label
        target = rows$feature,                  # taxon name
        logFC  = rows$logFC,
        qval   = rows$p_adjust,
        stringsAsFactors = FALSE
      )
    })

    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) stop("No edges produced after grouping by functions.")
    sankey_df <- do.call(rbind, pieces)

  } else {
    # Direct function -> taxon edges
    pieces <- lapply(names(pathway_sdaa_result), function(fid) {
      rows <- build_rows_from_df(pathway_sdaa_result[[fid]])
      if (is.null(rows)) return(NULL)
      data.frame(
        source = rep(fid, nrow(rows)),
        target = as.character(rows$feature),
        logFC  = rows$logFC,
        qval   = rows$p_adjust,
        stringsAsFactors = FALSE
      )
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) stop("No edges produced in function to taxon mode.")
    sankey_df <- do.call(rbind, pieces)
  }

  if (!nrow(sankey_df)) stop("No edges to plot after assembling sankey data.")

  # ---- optional source filter ---------------------------------------------------
  if (!is.null(taxa_list)) {
    keep_sources <- unique(as.character(taxa_list))
    mask <- tolower(sankey_df$source) %in% tolower(keep_sources)
    if (!any(mask)) stop("None of the requested sources in 'taxa_list' are present.")
    sankey_df <- sankey_df[mask, , drop = FALSE]
    if (!nrow(sankey_df)) stop("No edges remain after filtering by 'taxa_list'.")
  }

  # ---- nodes, links, colors -----------------------------------------------------
  nodes <- data.frame(name = unique(c(sankey_df$source, sankey_df$target)),
                      stringsAsFactors = FALSE)
  src_names <- unique(sankey_df$source)
  nodes$side <- ifelse(nodes$name %in% src_names, "left", "right")

  # Node/link groups for color scale
  nodes$node_group <- "All"
  sankey_df$link_group <- ifelse(sankey_df$logFC > 0, "Up", "Down")

  # Indices for networkD3 (0-based)
  sankey_df$source_id <- match(sankey_df$source, nodes$name) - 1L
  sankey_df$target_id <- match(sankey_df$target, nodes$name) - 1L

  # Link value: preserve your original exponential scaling
  sankey_df$value <- exp(sankey_df$logFC) * 10

  # Colour scale for nodes/links
  domain <- c("All","Up","Down")
  range  <- c(node_neutral, unname(node_color[c("Up","Down")]))

  color_scale <- htmlwidgets::JS(sprintf(
    "d3.scaleOrdinal().domain(%s).range(%s)",
    jsonlite::toJSON(domain, auto_unbox = TRUE),
    jsonlite::toJSON(range,  auto_unbox = TRUE)
  ))

  p <- networkD3::sankeyNetwork(
    Links       = sankey_df,
    Nodes       = nodes,
    Source      = "source_id",
    Target      = "target_id",
    Value       = "value",
    NodeID      = "name",
    LinkGroup   = "link_group",
    NodeGroup   = "node_group",
    colourScale = color_scale,
    fontSize    = 18,
    nodeWidth   = 30,
    sinksRight  = TRUE,
    margin      = margin,
    height      = height,
    width       = width
  )

  # ---- optional left-label wrapping/anchoring ----------------------------------
  if (isTRUE(label_left_sources)) {
    p <- htmlwidgets::onRender(p, "
      function(el, x){
        var fontSize = 18,  // keep in sync with R
            maxWidth = 220; // wrap width in px for LEFT labels

        var svg = d3.select(el).select('svg');

        // anchor & x-offsets
        svg.selectAll('.node text')
          .attr('text-anchor', function(d){
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? 'end' : 'start';
          })
          .attr('x', function(d){
            var w = (d.x1!==undefined && d.x0!==undefined) ? (d.x1 - d.x0)
                    : (d.dx!==undefined ? d.dx : 24);
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? -10 : (w + 10);
          });

        // wrap LEFT labels only
        svg.selectAll('.node')
          .filter(function(d){ return (d.side === 'left') || (d.x0 === 0 || d.x === 0); })
          .select('text')
          .each(function(d){
            var text = d3.select(this),
                words = (text.text() || '').split(/\\s+/).filter(Boolean),
                line = [], lineNumber = 0,
                lineHeight = 1.2, // ems
                y = (d.y0 + d.y1) / 2,
                x = +text.attr('x'),
                anchor = text.attr('text-anchor');

            text.text(null);
            var tspan = text.append('tspan')
                            .attr('x', x)
                            .attr('y', y)
                            .attr('dy', 0 + 'em')
                            .attr('text-anchor', anchor);

            for (var i=0; i<words.length; ++i) {
              line.push(words[i]);
              tspan.text(line.join(' '));
              if (tspan.node().getComputedTextLength() > maxWidth) {
                line.pop();
                tspan.text(line.join(' '));
                line = [words[i]];
                ++lineNumber;
                tspan = text.append('tspan')
                            .attr('x', x)
                            .attr('y', y)
                            .attr('dy', lineNumber * lineHeight + 'em')
                            .attr('text-anchor', anchor)
                            .text(words[i]);
              }
            }
          });

        // tooltip with full label
        svg.selectAll('.node').append('title').text(function(d){ return d.name; });
      }
    ")
  }

  p
}
