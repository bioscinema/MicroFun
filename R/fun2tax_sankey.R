#' Sankey diagram: functions (KO/EC/MetaCyc) → taxa/features
#'
#' Build an interactive Sankey diagram linking functional features (KO, EC, or
#' MetaCyc IDs) — optionally grouped to higher-level categories — to taxa or
#' features reported in differential analysis results (e.g., from
#' \code{pathway_sdaa()}).
#'
#' @description
#' The input \code{pathway_sdaa_result} must be a **named list** where names are
#' function IDs (KO/EC/MetaCyc), and each element is a data frame containing at
#' least \code{feature} (target/taxon), \code{logFC}, and one multiple-testing
#' column (e.g., \code{p_adjust}, \code{padj}, \code{qval}, \code{FDR}). The
#' function:
#' \enumerate{
#'   \item normalizes function IDs (e.g., strips \code{"ko:"}, \code{"ec:"}),
#'   \item (optionally) maps IDs to higher-level groups using reference tables
#'         shipped in \pkg{MicroFun} (\code{inst/extdata/*_reference.RData}),
#'   \item assembles edges \emph{source → target} (group or function → taxon),
#'   \item scales link width by \eqn{\exp(\mathrm{logFC}) \times 10},
#'   \item colors links by direction (Up/Down), and
#'   \item returns an interactive \pkg{networkD3} sankey widget.
#' }
#'
#' @param pathway_sdaa_result Named list. Names are function IDs (KO/EC/MetaCyc).
#'   Each element is a data frame with columns including:
#'   \itemize{
#'     \item \code{feature}: target/taxon (e.g., genus) label,
#'     \item \code{logFC}: log fold-change,
#'     \item a q/p-adjust column (any of: \code{p_adjust}, \code{padj},
#'           \code{qval}, \code{q_value}, \code{qvalue}, \code{p.adj}, \code{FDR}, \code{p_adj}).
#'   }
#' @param pathway Character. Functional namespace of the names in
#'   \code{pathway_sdaa_result}; one of \code{"KO"}, \code{"EC"}, \code{"MetaCyc"}.
#' @param group_functions Logical. If \code{TRUE}, map function IDs to higher
#'   groups before linking to taxa (KO → \code{PathwayL1}; EC/MetaCyc → \code{description}).
#'   Requires reference files in \pkg{MicroFun}'s \code{inst/extdata}.
#' @param taxa_list Optional character vector. If provided, keep only edges whose
#'   \emph{source} (group label or function ID) matches one of these entries.
#' @param node_color Named character vector of length 2 with colors for link
#'   groups: names must be \code{c("Up","Down")}. Default: orange for Up
#'   (\code{"#ffa551"}), blue for Down (\code{"#70afdf"}).
#' @param node_neutral Character. Color for nodes (neutral "All" category).
#'   Default \code{"#999999"}.
#' @param label_left_sources Logical. If \code{TRUE}, left-side node labels are
#'   right-aligned and wrapped via JavaScript for readability. Default \code{TRUE}.
#' @param margin Named list with pixel margins \code{top}, \code{right},
#'   \code{bottom}, \code{left}. Passed to \code{networkD3::sankeyNetwork()}.
#' @param height,width Numeric. Output size in pixels for the widget canvas.
#'
#'
#' @examples
#' \dontrun{
#' # Minimal toy example
#' sdaa_list <- list(
#'   K00001 = data.frame(feature = c("GenusA","GenusB"),
#'                       logFC = c( 1.1, -0.4),
#'                       p_adjust = c(0.02, 0.3)),
#'   K00002 = data.frame(feature = c("GenusB","GenusC"),
#'                       logFC = c(-1.6,  0.5),
#'                       p_adjust = c(0.01, 0.7))
#' )
#'
#' # Function → taxon (no grouping)
#' p1 <- fun2tax_sankey(sdaa_list, pathway = "KO", group_functions = FALSE)
#' p1
#'
#' # Group functions to higher-level categories (requires reference files)
#' p2 <- fun2tax_sankey(sdaa_list, pathway = "KO", group_functions = TRUE)
#' p2
#' }
#'
#'
#' @importFrom stats setNames
#' @importFrom jsonlite toJSON
#' @importFrom htmlwidgets JS onRender
#' @importFrom networkD3 sankeyNetwork
#' @export
fun2tax_sankey <- function(
    pathway_sdaa_result,                       # named list: names = function IDs; each df has feature, logFC, p_adjust
    pathway            = c("KO","EC","MetaCyc"),
    group_functions    = TRUE,                # if TRUE, map function IDs to higher-level groups (KO: PathwayL1; EC/MetaCyc: description)
    taxa_list          = NULL,                 # optional filter: keep only edges whose SOURCE (group or function) matches
    node_color         = c(Up = "#ffa551", Down = "#70afdf"),
    node_neutral       = "#999999",
    label_left_sources = TRUE,
    margin             = list(top=10, right=200, bottom=10, left=180),
    height             = 1200,
    width              = 1000
) {
  stopifnot(is.list(pathway_sdaa_result))
  pathway <- match.arg(pathway)

  # ---------- helpers ----------
  normalize_fun_ids <- function(ids, pathway = c("KO","EC","MetaCyc")) {
    pathway <- match.arg(pathway)
    ids <- trimws(as.character(ids))
    if (pathway == "KO") {
      ids <- toupper(gsub("^ko:", "", ids, ignore.case = TRUE))
    } else if (pathway == "EC") {
      ids <- gsub("^ec:", "", ids, ignore.case = TRUE)
      ids <- gsub("^EC[\\s_]*", "", ids, ignore.case = TRUE)
    } else {
      ids <- gsub("^metacyc:", "", ids, ignore.case = TRUE)
    }
    ids
  }

  load_reference_map <- function(pathway) {
    file_name <- switch(pathway,
                        KO      = "KO_reference.RData",
                        EC      = "EC_reference.RData",
                        MetaCyc = "MetaCyc_reference.RData"
    )
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
      setNames(as.character(ref$PathwayL1), as.character(normalize_fun_ids(ref$KO, "KO")))
    } else if (pathway == "EC") {
      setNames(as.character(ref$description), as.character(normalize_fun_ids(ref$EC, "EC")))
    } else {
      id_col <- if ("MetaCyc" %in% names(ref)) "MetaCyc" else if ("ID" %in% names(ref)) "ID" else NA
      if (is.na(id_col)) stop("MetaCyc reference must have an ID column")
      ids_norm <- normalize_fun_ids(ref[[id_col]], "MetaCyc")
      setNames(as.character(ref$description), as.character(ids_norm))
    }
  }

  get_qcol <- function(df) {
    cand <- c("p_adjust","padj","qval","q_value","qvalue","p.adj","FDR","p_adj")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  # ---------- normalize function IDs on the list ----------
  fun_ids_raw <- names(pathway_sdaa_result)
  if (is.null(fun_ids_raw)) stop("pathway_sdaa_result must be a *named* list (names = function IDs).")

  fun_ids_norm <- normalize_fun_ids(fun_ids_raw, pathway)
  keep <- !is.na(fun_ids_norm) & nzchar(fun_ids_norm)
  if (!any(keep)) stop("After normalization, no valid function IDs remain.")
  pathway_sdaa_result <- pathway_sdaa_result[keep]
  fun_ids_norm <- fun_ids_norm[keep]
  names(pathway_sdaa_result) <- fun_ids_norm

  # (Optional) load id -> group mapping
  id2group <- NULL
  if (isTRUE(group_functions)) {
    id2group <- load_reference_map(pathway)
  }

  # ---------- build sankey_df (function → taxon, or group → taxon) ----------
  if (isTRUE(group_functions)) {
    # group by higher-level label first
    grp_lab <- unname(id2group[fun_ids_norm])
    valid_fun <- !is.na(grp_lab) & nzchar(grp_lab)
    fun_ids_norm <- fun_ids_norm[valid_fun]
    grp_lab <- grp_lab[valid_fun]
    if (!length(fun_ids_norm)) stop("No function IDs matched the reference mapping after normalization.")
    groups <- split(fun_ids_norm, grp_lab)

    pieces <- lapply(names(groups), function(glab) {
      fids <- groups[[glab]]
      # collect rows from all member functions
      rows_list <- lapply(fids, function(fid) {
        df <- pathway_sdaa_result[[fid]]
        if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
        qcol <- get_qcol(df)
        if (is.null(qcol) || !"feature" %in% names(df) || !"logFC" %in% names(df)) return(NULL)
        keep <- !is.na(df$feature) & nzchar(df$feature) & !is.na(df[[qcol]])
        if (!any(keep)) return(NULL)
        data.frame(
          feature  = as.character(df$feature[keep]),
          logFC    = df$logFC[keep],
          p_adjust = df[[qcol]][keep],
          stringsAsFactors = FALSE
        )
      })
      rows_list <- Filter(Negate(is.null), rows_list)
      if (!length(rows_list)) return(NULL)
      rows <- do.call(rbind, rows_list)
      if (!nrow(rows)) return(NULL)

      # dedup per taxon within group: keep the row with the LARGER p_adjust (your rule)
      rows <- rows[order(rows$feature, rows$p_adjust, decreasing = TRUE), , drop = FALSE]
      rows <- rows[!duplicated(rows$feature), , drop = FALSE]

      data.frame(
        source = rep(glab, nrow(rows)),        # group label
        target = rows$feature,
        logFC  = rows$logFC,
        qval   = rows$p_adjust,
        stringsAsFactors = FALSE
      )
    })

    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) stop("No edges produced after grouping by functions.")
    sankey_df <- do.call(rbind, pieces)

  } else {
    # no grouping: function → taxon edges directly
    pieces <- lapply(names(pathway_sdaa_result), function(fid) {
      df <- pathway_sdaa_result[[fid]]
      if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
      qcol <- get_qcol(df)
      if (is.null(qcol) || !"feature" %in% names(df) || !"logFC" %in% names(df)) return(NULL)
      keep <- !is.na(df$feature) & nzchar(df$feature) & !is.na(df[[qcol]])
      if (!any(keep)) return(NULL)
      data.frame(
        source = rep(fid, sum(keep)),
        target = as.character(df$feature[keep]),
        logFC  = df$logFC[keep],
        qval   = df[[qcol]][keep],
        stringsAsFactors = FALSE
      )
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) stop("No edges produced in function→taxon mode.")
    sankey_df <- do.call(rbind, pieces)
  }

  if (nrow(sankey_df) == 0) stop("No edges to plot after assembling sankey data.")

  # ---------- optional source filter ----------
  if (!is.null(taxa_list)) {
    keep_sources <- unique(as.character(taxa_list))
    mask <- tolower(sankey_df$source) %in% tolower(keep_sources)
    if (!any(mask)) stop("None of the requested sources in 'taxa_list' are present.")
    sankey_df <- sankey_df[mask, , drop = FALSE]
  }

  # ---------- nodes, links, colors ----------
  nodes <- data.frame(
    name = unique(c(sankey_df$source, sankey_df$target)),
    stringsAsFactors = FALSE
  )
  src_names <- unique(sankey_df$source)
  nodes$side <- ifelse(nodes$name %in% src_names, "left", "right")

  # Separate node vs link groups so neutral node color is effective
  nodes$node_group <- "All"
  sankey_df$link_group <- ifelse(sankey_df$logFC > 0, "Up", "Down")

  sankey_df$source_id <- match(sankey_df$source, nodes$name) - 1L
  sankey_df$target_id <- match(sankey_df$target, nodes$name) - 1L
  sankey_df$value     <- exp(sankey_df$logFC)*10

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
    height = height, width = width
  )

  if (label_left_sources) {
    p <- htmlwidgets::onRender(p, "
    function(el, x){
      var fontSize = 18,  // keep in sync with R
          maxWidth = 220; // wrap width in px for LEFT labels

      var svg = d3.select(el).select('svg');

      // set anchors & x offsets (as you had)
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

      // tooltip with full label on hover
      svg.selectAll('.node').append('title').text(function(d){ return d.name; });
    }
  ")
  }


  p
}
