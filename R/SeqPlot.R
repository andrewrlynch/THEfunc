#' Load cytoband table for ideogram rendering
#'
#' This function loads a UCSC-style cytoband file for ideogram drawing.
#' If no path is given, it loads the built-in cytoband_hg38 dataset from the package.
#'
#' @param path Optional path to a UCSC-style cytoband table (TSV format).
#'
#' @return A data.frame with cytoband annotations.
#' @export
loadCytobands <- function(path = NULL) {
  if (is.null(path)) {
    data("cytoband_hg38", package = "THEfunc", envir = environment())
    return(cytoband_hg38)
  } else {
    return(read.table(
      path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    ))
  }
}



#' defaultGenomeWindows
#' @name defaultGenomeWindows
#' @description
#' Generates a GRanges object of all standard hg38 chromosomes to serve as the default plotting windows.
#' @export
#' @author Andrew Lynch
defaultGenomeWindows <- function() {
  chr_lengths <- c(
    "chr1"  = 248956422,
    "chr2"  = 242193529,
    "chr3"  = 198295559,
    "chr4"  = 190214555,
    "chr5"  = 181538259,
    "chr6"  = 170805979,
    "chr7"  = 159345973,
    "chr8"  = 145138636,
    "chr9"  = 138394717,
    "chr10" = 133797422,
    "chr11" = 135086622,
    "chr12" = 133275309,
    "chr13" = 114364328,
    "chr14" = 107043718,
    "chr15" = 101991189,
    "chr16" = 90338345,
    "chr17" = 83257441,
    "chr18" = 80373285,
    "chr19" = 58617616,
    "chr20" = 64444167,
    "chr21" = 46709983,
    "chr22" = 50818468,
    "chrX"  = 156040895,
    "chrY"  = 57227415
  )

  GRanges(seqnames = names(chr_lengths),
          ranges = IRanges(start = 1, end = chr_lengths))
}



#' Transform Genomic X Coordinates to Absolute Plot X Coordinates
#'
#' Converts genomic x coordinates into absolute plot x coordinates using the layout information.
#'
#' @param x A numeric vector of genomic positions.
#' @param w An integer vector of window indices corresponding to each genomic position in x.
#' @param layout A list returned by setupGlobalLayout() that includes global_window_starts, global_window_ends,
#'   and global_windows.
#'
#' @return A numeric vector of absolute x coordinates.
#' @export
globalTransformX <- function(x, w, layout) {

  win <- layout$global_windows[w]
  win_start <- start(win)
  win_width <- width(win)

  x_rel <- (x - win_start) / win_width
  x_rel <- pmax(pmin(x_rel, 1), 0)

  abs_x <- unname(layout$global_window_starts)[w] +
    x_rel * (unname(layout$global_window_ends)[w] - unname(layout$global_window_starts)[w])
  return(abs_x)
}




#' Transform Relative Y Coordinates to Absolute Plot Y Coordinates
#'
#' Converts relative y coordinates (values typically between 0 and 1) into absolute plot y coordinates
#' using the layout information. If per-track y ranges (track_y_min / track_y_max) are provided, they are used for scaling.
#'
#' @param y A numeric vector of relative y values.
#' @param t An integer vector of track indices corresponding to each y value.
#' @param layout A list returned by setupGlobalLayout() that includes global_track_bottoms and global_track_tops.
#' @param track_y_min Optional numeric vector of minimum y values for each track. Defaults to 0 for all tracks.
#' @param track_y_max Optional numeric vector of maximum y values for each track. Defaults to 1 for all tracks.
#'
#' @return A numeric vector of absolute y coordinates.
#' @export
globalTransformY <- function(y, t, layout_info, track_y_min, track_y_max) {
  if (length(y) == 0 || length(t) == 0) return(numeric(0))

  stopifnot(length(y) == length(t))  # Sanity check

  rel <- numeric(length(y))          # Initialize relative y
  abs_y <- numeric(length(y))        # Final absolute y

  for (i in seq_along(y)) {
    track_idx <- t[i]
    # Fallback if track index is missing or out of bounds
    if (is.na(track_idx) || track_idx < 1 || track_idx > length(layout_info$global_track_bottoms)) {
      abs_y[i] <- NA_real_
      next
    }

    ymin <- track_y_min[track_idx]
    ymax <- track_y_max[track_idx]

    if (ymin == ymax) {
      rel[i] <- 0.5
    } else {
      rel[i] <- (y[i] - ymin) / (ymax - ymin)
    }

    rel[i] <- max(0, min(1, rel[i]))  # Clamp to [0,1]

    bottom <- layout_info$global_track_bottoms[track_idx]
    top <- layout_info$global_track_tops[track_idx]
    abs_y[i] <- bottom + rel[i] * (top - bottom)
  }

  return(abs_y)
}




#' Set Up Global Layout for SeqPlot
#'
#' Computes the plotting canvas layout including track positions, window positions,
#' and scaling factors based on the SeqPlot parameters.
#'
#' @param sp A SeqPlot object containing tracks, track gaps, track heights,
#'   globalWindows, windowGap, and margin.
#' @param plotWidth Numeric. The width of the plotting area (excluding margins).
#'   Default is 1000.
#' @param globalWindows Optional GRanges object defining custom genomic windows.
#'
#' @return A list containing:
#'   \item{margin}{The margin list from sp.}
#'   \item{global_track_heights}{Numeric vector of absolute pixel heights for each track.}
#'   \item{global_track_bottoms}{Numeric vector of the bottom y-positions of each track.}
#'   \item{global_track_tops}{Numeric vector of the top y-positions of each track.}
#'   \item{global_plot_height}{Total height of the plotting canvas including margins.}
#'   \item{global_window_starts}{Numeric vector of the absolute x-coordinate start for each window.}
#'   \item{global_window_ends}{Numeric vector of the absolute x-coordinate end for each window.}
#'   \item{scaleFactors}{Numeric vector of scaling factors for each global window.}
#'   \item{global_windows}{The GRanges object of global windows.}
#'
#' @export
setupGlobalLayout <- function(sp, plotWidth = 1000, globalWindows = NULL) {
  margin <- sp@margin
  totalTracks <- length(sp@tracks)

  if (is.null(globalWindows)) {
    global_windows <- defaultGenomeWindows()
  } else {
    global_windows <- globalWindows
  }

  # Order global_windows by seqnames and start
  windowOrder <- unique(as.character(seqnames(global_windows)))
  global_windows <- global_windows[
    as.character(seqnames(global_windows)) %in% windowOrder
  ]
  global_windows <- global_windows[
    order(factor(as.character(seqnames(global_windows)), levels = windowOrder), start(global_windows))
  ]
  num_windows <- length(global_windows)

  # Scale factors (default to Mb if not defined)
  if (is.null(mcols(global_windows)$scale)) {
    scale_factors <- rep(1e-6, num_windows)
  } else {
    scale_factors <- mcols(global_windows)$scale
    scale_factors[is.na(scale_factors)] <- 1e-6
  }

  # Genomic widths scaled
  effective_widths <- width(global_windows) * scale_factors
  total_effective_width <- sum(effective_widths)

  assigned_widths <- (effective_widths / total_effective_width) *
    (plotWidth - (num_windows - 1) * sp@windowGap)

  global_x_start <- numeric(num_windows)
  global_x_end <- numeric(num_windows)
  global_x_start[1] <- margin$left
  global_x_end[1] <- global_x_start[1] + assigned_widths[1]
  if (num_windows > 1) {
    for (i in 2:num_windows) {
      global_x_start[i] <- global_x_end[i - 1] + sp@windowGap
      global_x_end[i] <- global_x_start[i] + assigned_widths[i]
    }
  }

  # Track heights (proportional)
  if (is.null(sp@trackHeights)) {
    sp@trackHeights <- rep(1, totalTracks)
  }
  relative_heights <- sp@trackHeights / sum(sp@trackHeights)
  available_height <- 1000 - (totalTracks - 1) * sp@trackGap
  track_pixel_heights <- relative_heights * available_height

  global_track_bottom <- numeric(totalTracks)
  global_track_top <- numeric(totalTracks)
  global_track_bottom[1] <- margin$bottom
  global_track_top[1] <- global_track_bottom[1] + track_pixel_heights[1]
  if (totalTracks > 1) {
    for (i in 2:totalTracks) {
      global_track_bottom[i] <- global_track_top[i - 1] + sp@trackGap
      global_track_top[i] <- global_track_bottom[i] + track_pixel_heights[i]
    }
  }

  totalWidth <- plotWidth + margin$left + margin$right
  totalHeight <- global_track_top[totalTracks] + margin$top

  list(
    global_windows = global_windows,
    global_window_starts = global_x_start,
    global_window_ends = global_x_end,
    global_x_start = global_x_start,
    global_x_end = global_x_end,
    global_windows = global_windows,
    global_scale_factors = scale_factors,
    global_track_bottoms = global_track_bottom,
    global_track_tops = global_track_top,
    plot_origin_x = margin$left,
    plot_origin_y = margin$bottom,
    global_margins = margin,
    global_plot_width = totalWidth,
    global_plot_height = totalHeight,
    global_margin = margin
  )
}





# SeqLink drawing ----

#' bezier3
#' @name bezier3
#' @description
#' Function that returns coordinates of a cubic bezier curve
#' @export
#' @author Andrew Lynch
bezier3 <- function(P, t) {
  n <- 3
  B <- sapply(0:n, function(i) choose(n, i) * (1 - t)^(n - i) * t^i)
  x <- colSums(B * P[,1])
  y <- colSums(B * P[,2])
  cbind(x, y)
}



#' bezier4
#' @name bezier4
#' @description
#' Function that returns coordinates of a cubic bezier curve
#' @export
#' @author Andrew Lynch
bezier4 <- function(P, t) {
  n <- 4
  B <- sapply(0:n, function(i) choose(n, i) * (1 - t)^(n - i) * t^i)
  x <- colSums(B * P[,1])
  y <- colSums(B * P[,2])
  cbind(x, y)
}



#' bezier5
#' @name bezier5
#' @description
#' Function that returns coordinates of a 5th order bezier curve
#' @export
#' @author Andrew Lynch
bezier5 <- function(P, t) {
  n <- 5
  B <- sapply(0:n, function(i) choose(n, i) * (1 - t)^(n - i) * t^i)
  x <- colSums(B * P[,1])
  y <- colSums(B * P[,2])
  cbind(x, y)
}



#' bezier6
#' @name bezier6
#' @description
#' Function that returns coordinates of a 6th order bezier curve
#' @export
#' @author Andrew Lynch
bezier6 <- function(P, t) {
  n <- 6
  B <- sapply(0:n, function(i) choose(n, i) * (1 - t)^(n - i) * t^i)
  x <- colSums(B * P[, 1])
  y <- colSums(B * P[, 2])
  cbind(x, y)
}

# Seq Classes ----

#' SeqFeature-Class
#' @name SeqFeature-Class
#' @description
#' S4 class for SeqFeatures
#' @author Andrew Lynch
#' @exportClass SeqFeature
setClass("SeqFeature",
  slots = list(
    gr = "GRanges",       # Original genomic ranges
    x0 = "numeric",       # Mapped plot start (genomic coordinate)
    x1 = "numeric",       # Mapped plot end (genomic coordinate)
    y0 = "numeric",       # Mapped y coordinate (lower)
    y1 = "numeric",       # Mapped y coordinate (upper)
    type = "character",   # Type of feature (e.g., "point", "segment", etc.)
    color = "character"   # Color for plotting the feature
  )
)
#' SeqFeature
#'
#' This function constructs a SeqFeature for plotting genomic features.
#'
#' @param gr A GRanges object with genomic coordinates.
#' @param type Character string indicating the feature type: "segment", "point", or "box".
#' @param yCol Metadata column name for y0 values.
#' @param y1Col Metadata column name for y1 values.
#' @param colorCol Metadata column name for per-feature color (optional).
#' @param color Default color if no column is specified.
#' @return A SeqFeature object.
#' @export
SeqFeature <- function(gr, type = c("segment", "point", "bar"),
                       yCol = NULL, y1Col = NULL,
                       colorCol = NULL, color = "black", ...) {
  type <- match.arg(type)
  x0 <- start(gr)
  x1 <- if (type == "segment" || type == "bar") end(gr) else start(gr)

  if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
    y0 <- mcols(gr)[[yCol]]
  } else {
    y0 <- if (type == "bar") rep(0, length(gr)) else rep(NA_real_, length(gr))
  }

  if (!is.null(y1Col) && y1Col %in% names(mcols(gr))) {
    y1 <- mcols(gr)[[y1Col]]
  } else {
    y1 <- if (type == "bar") rep(1, length(gr)) else y0
  }

  # Colors
  if (!is.null(colorCol) && colorCol %in% names(mcols(gr))) {
    colors <- as.character(mcols(gr)[[colorCol]])
  } else {
    colors <- rep(color, length(gr))
  }

  new("SeqFeature",
      gr = gr,
      x0 = as.numeric(x0),
      x1 = as.numeric(x1),
      y0 = as.numeric(y0),
      y1 = as.numeric(y1),
      type = type,
      color = colors)
}

## SeqFeature Subclasses ----

#' @title SeqPoint Class
#' @description Subclass of SeqFeature for point features (e.g., SNPs)
#' @exportClass SeqPoint
setClass("SeqPoint", contains = "SeqFeature")

#' @title SeqSegment Class
#' @description Subclass of SeqFeature for segment features (e.g., genes)
#' @exportClass SeqSegment
setClass("SeqSegment", contains = "SeqFeature")

#' @title SeqBar Class
#' @description Subclass of SeqFeature for rectangular blocks (e.g., bars, cytobands)
#' @exportClass SeqBar
setClass("SeqBar", contains = "SeqFeature")

#' @title SeqPoint
#' @description Constructs a SeqPoint object for plotting point features.
#' @param gr GRanges object
#' @param yCol Metadata column for y-axis value
#' @param colorCol Metadata column for color (optional)
#' @param color Default color
#' @export
SeqPoint <- function(gr, yCol = NULL, colorCol = NULL, color = "black") {
  sf <- SeqFeature(gr, type = "point", yCol = yCol, colorCol = colorCol, color = color)
  new("SeqPoint", sf)
}

#' @title SeqSegment
#' @description Constructs a SeqSegment object for plotting genomic segments.
#' @param gr GRanges object
#' @param yCol Metadata column for start y-value
#' @param y1Col Metadata column for end y-value
#' @param colorCol Metadata column for color (optional)
#' @param color Default color
#' @export
SeqSegment <- function(gr, yCol = NULL, y1Col = NULL, colorCol = NULL, color = "black") {
  sf <- SeqFeature(gr, type = "segment", yCol = yCol, y1Col = y1Col, colorCol = colorCol, color = color)
  new("SeqSegment", sf)
}

#' @title SeqBar
#' @description Constructs a SeqBar object for plotting rectangular blocks.
#' @param gr GRanges object
#' @param yCol Metadata column for bottom y-value
#' @param y1Col Metadata column for top y-value
#' @param colorCol Metadata column for color (optional)
#' @param color Default color
#' @export
SeqBar <- function(gr, yCol = NULL, y1Col = NULL, colorCol = NULL, color = "black") {
  sf <- SeqFeature(gr, type = "bar", yCol = yCol, y1Col = y1Col, colorCol = colorCol, color = color)
  new("SeqBar", sf)
}




# SeqLink ----------------------------------------------------------------------

#' SeqLink-Class
#' @name SeqLink-Class
#' @description
#' S4 class for SeqLinks
#' @author Andrew Lynch
#' @exportClass SeqLink
setClass("SeqLink",
         slots = list(
           gr1 = "GRanges",       # Start positions
           gr2 = "GRanges",       # End positions
           color = "character",   # Colors per link
           y = "numeric",         # Optional height scale or absolute y-value
           y0 = "numeric",
           y1 = "numeric",
           t0 = "numeric",
           t1 = "numeric",
           orientation = "character"
         ))


#' SeqLink
#' @description Constructs a base SeqLink object.
#' @export
SeqLink <- function(gr1, gr2,
                    color = "black", colorCol = NULL,
                    y = NA_real_, yCol = NULL,
                    y0 = NA_real_, y0Col = NULL,
                    y1 = NA_real_, y1Col = NULL,
                    t0 = 0, t1 = 0,
                    orientation = "*") {

  if (length(gr1) != length(gr2)) stop("gr1 and gr2 must have same length.")
  n <- length(gr1)

  # Metadata parsing
  y_vals <- if (!is.null(yCol) && yCol %in% names(mcols(gr1))) mcols(gr1)[[yCol]] else rep(y, n)
  y0_vals <- if (!is.null(y0Col) && y0Col %in% names(mcols(gr1))) mcols(gr1)[[y0Col]] else rep(y0, n)
  y1_vals <- if (!is.null(y1Col) && y1Col %in% names(mcols(gr2))) mcols(gr2)[[y1Col]] else rep(y1, n)
  colors <- if (!is.null(colorCol) && colorCol %in% names(mcols(gr1))) as.character(mcols(gr1)[[colorCol]]) else rep(color, n)

  new("SeqLink",
      gr1 = gr1,
      gr2 = gr2,
      color = colors,
      y = as.numeric(y_vals),
      y0 = as.numeric(y0_vals),
      y1 = as.numeric(y1_vals),
      t0 = rep(t0, n),
      t1 = rep(t1, n),
      orientation = rep(orientation, n))
}


drawLinkStub <- function(x, y, direction = "right",
                         layout = NULL, track_idx = NULL,
                         y_top = NULL, height = NULL,
                         track_y_min = NULL, track_y_max = NULL,
                         fixedStubs = TRUE,
                         stub_len = 25, col = "gray30", lwd = 1.5) {

  # Determine y_top logic
  if (is.null(y_top)) {
    if (!is.null(layout) && !is.null(track_idx)) {
      if (fixedStubs) {
        y_top <- layout$global_track_tops[track_idx]
      } else if (!is.null(height)) {
        y_top <- globalTransformY(height, track_idx, layout,
                                  track_y_min = rep(track_y_min, length(layout$global_track_bottoms)),
                                  track_y_max = rep(track_y_max, length(layout$global_track_bottoms)))
      } else {
        stop("Need 'height' if fixedStubs = FALSE.")
      }
    } else if (!is.null(height)) {
      y_top <- y + height
    } else {
      stop("Could not determine y_top: provide layout, height, or y_top.")
    }
  }

  # 1. Vertical stem
  segments(x, y, x, y_top, col = col, lwd = lwd)

  # 2. Hook (45° curve)
  n <- 20
  t <- seq(0, 45 * pi / 180, length.out = n)
  if (direction == "right") {
    cx <- x + stub_len * (1 - cos(t))
    cy <- y_top + stub_len * sin(t)
  } else if (direction == "left") {
    cx <- x - stub_len * (1 - cos(t))
    cy <- y_top + stub_len * sin(t)
  }


  lines(cx, cy, col = col, lwd = lwd)
  #arrows(tail(cx, 1), tail(cy, 1), tail(cx, 1) + ifelse(direction == "right", 5, -5), tail(cy, 1), length = 0.05, col = col, lwd = lwd)
}






## SeqLink Subclasses ----------------------------------------------------------
### SeqArc ----------------------------------------------------------

#' SeqArc-Class
#' @description S4 subclass of SeqLink for curved arcs with optional height or orientation.
#' @exportClass SeqArc
setClass("SeqArc", contains = "SeqLink")


#' SeqArc
#' @description Constructs a SeqArc object for curved connections.
#' @export
SeqArc <- function(gr1, gr2,
                   color = "black", colorCol = NULL,
                   y = NA_real_, yCol = NULL,
                   y0 = NA_real_, y0Col = NULL,
                   y1 = NA_real_, y1Col = NULL,
                   t0 = 0, t1 = 0,
                   orientation = "*") {

  sf <- SeqLink(gr1 = gr1, gr2 = gr2,
                color = color, colorCol = colorCol,
                y = y, yCol = yCol,
                y0 = y0, y0Col = y0Col,
                y1 = y1, y1Col = y1Col,
                t0 = t0, t1 = t1,
                orientation = orientation)

  new("SeqArc", sf)
}





#' drawArc
#' @description Draws a cubic Bezier arc between two points.
#' @param x0, y0 Start coordinates
#' @param x1, y1 End coordinates
#' @param orientation Arc orientation: "*", "up", "down", "+", "-", etc.
#' @param col Color
#' @param lwd Line width
#' @export
drawArc <- function(x0, y0, x1, y1, orientation = "*", col = "black", lwd = 2) {
  upward <- orientation %in% c("*", "up", "+", "-/+", "-/-")
  sign <- ifelse(upward, 1, -1)

  # Midpoint
  mx <- (x0 + x1) / 2
  my <- (y0 + y1) / 2

  # Arc height is a fraction of span
  dx <- abs(x1 - x0)
  arc_offset <- sign * dx * 0.25

  # Control point placed perpendicular to line connecting start–end
  ctrl_x <- mx
  ctrl_y <- my + arc_offset

  # Cubic Bezier control points (start, ctrl, end)
  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^2 * x0 + 2 * (1 - t) * t * ctrl_x + t^2 * x1
  bez_y <- (1 - t)^2 * y0 + 2 * (1 - t) * t * ctrl_y + t^2 * y1

  lines(bez_x, bez_y, col = col, lwd = lwd)
}

### SeqArch ----------------------------------------------------------
#' SeqArch-Class
#' @description Bounded arch link with vertical lines and an arc on top.
#' @exportClass SeqArch
setClass("SeqArch",
         contains = "SeqLink",
         slots = list(
           height = "numeric",  # apex y-position (in track-relative units, e.g., 0–1)
           curve = "ANY"
         ))


#' SeqArch
#' @description Constructs a SeqArch object with vertical lines and an arc above.
#' @export
SeqArch <- function(gr1, gr2,
                    color = "black", colorCol = NULL,
                    y0 = 0, y0Col = NULL,
                    y1 = 0, y1Col = NULL,
                    t0 = 0, t1 = 0,
                    height = 1,
                    orientation = "*",
                    curve = "length") {

  sf <- SeqLink(gr1 = gr1, gr2 = gr2,
                color = color, colorCol = colorCol,
                y0 = y0, y0Col = y0Col,
                y1 = y1, y1Col = y1Col,
                t0 = t0, t1 = t1,
                orientation = orientation)

  n <- length(gr1)

  new("SeqArch", sf, height = rep(height, n), curve = rep(curve, n))
}



#' drawSeqArch
#' @description Draws a bounded arch with vertical lines and a curved top.
#' @param x0, y0 Start base (bottom of left stem)
#' @param x1, y1 End base (bottom of right stem)
#' @param top0, top1 Absolute Y positions of the top of the stems (in arc track)
#' @param orientation Character: "up", "down", "+", "-", "*", etc.
#' @param col Color
#' @param lwd Line width
#' @export
drawSeqArch <- function(x0, y0, x1, y1, top0, top1,
                        orientation = "*", curve = "length",
                        col = "black", lwd = 2) {

  upward <- orientation %in% c("*", "up", "+", "-/+", "-/-")
  sign <- ifelse(upward, 1, -1)

  span <- abs(x1 - x0)

  # Determine vertical arc "bulge" height
  if (is.numeric(curve)) {
    curve_offset <- sign * 100 * curve
  } else if (curve == "equal") {
    curve_offset <- sign * 100 * 0.25  # default uniform curve
  } else if (curve == "length") {
    curve_offset <- sign * span * 0.25  # scale with span (same multiplier)
  } else {
    warning("Unknown curve value in drawSeqArch(); using default 0.25")
    curve_offset <- sign * span * 0.25
  }

  # Controls how "broad" or flat the arch appears (adjustable factor)
  ctrl_spread <- 0.1  # increase for broader arches (e.g., 0.35, 0.5)

  # Horizontal offset from endpoints
  dx <- abs(x1 - x0)
  ctrl_dx <- dx * ctrl_spread

  mid_y <- max(top0, top1)
  ctrl_x1 <- x0 + ctrl_dx
  ctrl_x2 <- x1 - ctrl_dx

  ctrl_y1 <- mid_y + curve_offset
  ctrl_y2 <- mid_y + curve_offset

  P0 <- c(x0, top0)
  P1 <- c(ctrl_x1, ctrl_y1)
  P2 <- c(ctrl_x2, ctrl_y2)
  P3 <- c(x1, top1)

  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^3 * P0[1] + 3 * (1 - t)^2 * t * P1[1] +
    3 * (1 - t) * t^2 * P2[1] + t^3 * P3[1]
  bez_y <- (1 - t)^3 * P0[2] + 3 * (1 - t)^2 * t * P1[2] +
    3 * (1 - t) * t^2 * P2[2] + t^3 * P3[2]

  lines(bez_x, bez_y, col = col, lwd = lwd)

  # Draw vertical stems
  segments(x0, y0, x0, top0, col = col, lwd = lwd)
  segments(x1, y1, x1, top1, col = col, lwd = lwd)


}



### SeqRecon ----
### Inspired by GH: cortes-ciriano-lab/ReConPlot


#' SeqRecon-Class
#' @description ReCon-style rearrangement arcs, subclass of SeqArch with automatic tiering and coloring.
#' @exportClass SeqRecon
setClass("SeqRecon", contains = "SeqArch")

#' SeqRecon
#' @description Constructs a ReCon-style arc from paired breakpoints.
#' @export
SeqRecon <- function(gr1, gr2,
                     orientation = "*", orientationCol = NULL,
                     color = NULL, colorCol = NULL,
                     y0 = 0, y1 = 0, t0 = 0, t1 = 0,
                     curve = "equal") {

  n <- length(gr1)

  # Extract orientation from metadata if needed
  if (!is.null(orientationCol) && orientationCol %in% names(mcols(gr1))) {
    orientation_vals <- mcols(gr1)[[orientationCol]]
  } else if (length(orientation) == 1) {
    orientation_vals <- rep(orientation, n)
  } else {
    orientation_vals <- orientation  # already a vector of correct length
  }

  # Infer SV type
  sv_type <- mapply(function(ori, chr1, chr2) {
    if (chr1 != chr2) return("TRA")
    switch(ori,
           "+/+" = "h2hINV", "-/-" = "t2tINV",
           "+/-" = "DEL", "-/+" = "DUP",
           "unknown")
  }, orientation_vals, as.character(seqnames(gr1)), as.character(seqnames(gr2)))

  # Define height per type (these are Y positions in the arc track)
  sv_height <- c(h2hINV = 0.25, t2tINV = 0.25, DEL = 0.5, DUP = 0.5, TRA = 0.75)
  heights <- sv_height[sv_type]
  heights[is.na(heights)] <- 0.5

  # Define default colors
  sv_color <- c(h2hINV = "#FF7800", t2tINV = "#FFB900",
                DEL = "#0E97CF", DUP = "#E92F2F",
                TRA = "#74B21A")
  colors <- if (!is.null(colorCol) && colorCol %in% names(mcols(gr1))) {
    as.character(mcols(gr1)[[colorCol]])
  } else {
    sv_color[sv_type]
  }

  new("SeqRecon",
      SeqArch(gr1 = gr1, gr2 = gr2,
              orientation = orientation_vals,
              y0 = y0, y1 = y1,
              t0 = t0, t1 = t1,
              height = heights,
              curve = curve,
              color = colors))
}



### SeqString ----
### Inspired by GH: mskilab-org/gTrack and mskilab-org/gGnome

#' SeqString-Class
#' @description A smooth, signed Bezier curve between two genomic points with horizontal entry/exit
#' @exportClass SeqString
setClass("SeqString", contains = "SeqLink")

#' SeqString
#' @description Constructs a SeqString object from a BED-like input
#' @param df A data.frame or DataFrame with columns chr0, x0, y0, chr1, x1, y1, orientation
#' @export
SeqString <- function(df,
                      t0 = 0, t1 = 0,
                      color = "black", colorCol = NULL,
                      orientation = "*", orientationCol = "orientation") {

  if (!all(c("chr0", "x0", "y0", "chr1", "x1", "y1") %in% names(df))) {
    stop("Input must contain chr0, x0, y0, chr1, x1, y1 columns.")
  }

  n <- nrow(df)

  # Orientation from column if available
  if (!is.null(orientationCol) && orientationCol %in% names(df)) {
    ori <- as.character(df[[orientationCol]])
  } else if (length(orientation) == 1) {
    ori <- rep(orientation, n)
  } else {
    ori <- orientation
  }

  # Convert to GRanges (we'll use gr1 and gr2 as placeholders)
  gr1 <- GRanges(seqnames = df$chr0,
                 ranges = IRanges(start = df$x0, width = 1),
                 y0 = df$y0)

  gr2 <- GRanges(seqnames = df$chr1,
                 ranges = IRanges(start = df$x1, width = 1),
                 y1 = df$y1)

  # Colors
  col_vals <- if (!is.null(colorCol) && colorCol %in% names(df)) {
    as.character(df[[colorCol]])
  } else {
    rep(color, n)
  }

  new("SeqString", SeqLink(gr1 = gr1, gr2 = gr2,
                           t0 = rep(t0, n), t1 = rep(t1, n),
                           y0 = df$y0, y1 = df$y1,
                           orientation = ori,
                           color = col_vals))
}



#' drawSeqString
#' @description Draws a 6-point cubic Bézier DNA string curve
#' @export
drawSeqString <- function(x0, y0, x1, y1, orientation = "*",
                          col = "black", lwd = 1.5, buffer = NULL) {

  dx <- abs(x1 - x0)
  dy <- abs(y1 - y0)

  # Set horizontal "stretch" buffer
  if (is.null(buffer)) {
    buffer <- dx * 0.15
  }

  # Control points based on orientation
  switch(orientation,
         "+/+" = {
           P0 <- c(x0, y0)
           P1 <- c(x0 + buffer, y0)
           P4 <- c(x1 - buffer, y1)
           P5 <- c(x1, y1)
         },
         "-/-" = {
           P0 <- c(x0, y0)
           P1 <- c(x0 - buffer, y0)
           P4 <- c(x1 + buffer, y1)
           P5 <- c(x1, y1)
         },
         "+/-" = {
           P0 <- c(x0, y0)
           P1 <- c(x0 + buffer, y0)
           P4 <- c(x1 + buffer, y1)
           P5 <- c(x1, y1)
         },
         "-/+" = {
           P0 <- c(x0, y0)
           P1 <- c(x0 - buffer, y0)
           P4 <- c(x1 - buffer, y1)
           P5 <- c(x1, y1)
         },
         "*" = {
           # Default upward arc (like SeqArch)
           ctrl_x <- (x0 + x1) / 2
           ctrl_y <- max(y0, y1) + dx * 0.25
           t <- seq(0, 1, length.out = 100)
           bez_x <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * ctrl_x +
             3 * (1 - t) * t^2 * ctrl_x + t^3 * x1
           bez_y <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * ctrl_y +
             3 * (1 - t) * t^2 * ctrl_y + t^3 * y1
           lines(bez_x, bez_y, col = col, lwd = lwd)
           return(invisible())
         },
         {
           warning(sprintf("Unknown orientation '%s'; skipping draw", orientation))
           return(invisible())
         }
  )

  # Intermediate control points for smoothness
  P2 <- P1 + (P4 - P1) * 0.4
  P3 <- P1 + (P4 - P1) * 0.6

  # Generate Bezier curve
  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^5 * P0[1] +
    5 * (1 - t)^4 * t * P1[1] +
    10 * (1 - t)^3 * t^2 * P2[1] +
    10 * (1 - t)^2 * t^3 * P3[1] +
    5 * (1 - t) * t^4 * P4[1] +
    t^5 * P5[1]

  bez_y <- (1 - t)^5 * P0[2] +
    5 * (1 - t)^4 * t * P1[2] +
    10 * (1 - t)^3 * t^2 * P2[2] +
    10 * (1 - t)^2 * t^3 * P3[2] +
    5 * (1 - t) * t^4 * P4[2] +
    t^5 * P5[2]

  lines(bez_x, bez_y, col = col, lwd = lwd)
}



# SeqAnnotation ----

#' SeqAnnotation-Class
#'
#' Generic annotation class for drawing composite genomic annotations (e.g., ideograms, gene models).
#' @slot annotationType Type of annotation (e.g., "ideogram").
#' @slot data Internal data required for rendering the annotation.
#' @exportClass SeqAnnotation
setClass("SeqAnnotation",
         slots = list(
           annotationType = "character",
           y0 = "numeric",
           y1 = "numeric"
         ))




#' SeqAnnotation Constructor
#'
#' Create a SeqAnnotation object.
#' @export
setClass("SeqAnnotation",
         slots = list(annotationType = "character"))


## AnnoFeature ----
#' AnnoFeature-Class
#' @description Basic annotation box with optional labels and shapes.
#' @export
setClass("AnnoFeature",
         contains = "SeqAnnotation",
         slots = list(
           gr = "GRanges",
           labels = "character",
           color = "character",
           orientation = "character",
           shape = "character",
           x0 = "numeric",
           x1 = "numeric",
           y0 = "numeric",
           y1 = "numeric"
         )
)





#' AnnoFeature
#' @description Constructs an AnnoFeature object from gene or enhancer annotations
#' @param gr A GRanges object with annotation features (e.g., from RefSeq)
#' @param labelCol Column name for feature labels (e.g., gene names)
#' @param colorCol Optional column for color by feature type or category
#' @export
AnnoFeature <- function(gr,
                        labelCol = NULL,
                        color = "gray30",
                        colorCol = NULL,
                        orientation = "*",
                        shape = "rect",
                        x0 = start(gr),
                        x1 = end(gr),
                        y0 = 0,
                        y1 = 1) {

  if (!inherits(gr, "GRanges")) stop("Input must be a GRanges object.")

  labels <- if (!is.null(labelCol) && labelCol %in% names(mcols(gr))) {
    as.character(mcols(gr)[[labelCol]])
  } else {
    rep(NA_character_, length(gr))
  }

  colors <- if (!is.null(colorCol) && colorCol %in% names(mcols(gr))) {
    as.character(mcols(gr)[[colorCol]])
  } else {
    rep(color, length(gr))
  }
  n <- length(gr)
  shapes <- rep(shape, n)
  y0 = rep(y0, n)
  y1 = rep(y1,n)
  strand_vals <- as.character(strand(gr))
  orientations <- ifelse(strand_vals %in% c("+", "-"), strand_vals, "*")

  new("AnnoFeature", gr = gr, labels = labels, color = colors, shape = shapes, orientation = orientations, y0 = y0, y1 = y1, x0 = as.numeric(x0), x1 = as.numeric(x1))
}


#' drawAnnoFeature
#' @description Draw annotation features with different shape styles
#' @export
drawAnnoFeature <- function(gr, labels = NULL, orientation = NULL,
                            shape = "rect",
                            x0 = NULL, x1 = NULL,
                            xlim = NULL, ybase = 0, tier_height = 1,
                            y0s = NULL, y1s = NULL,
                            col = "gray30", border = NA, label_cex = 0.6,
                            label_offset = 4,
                            returnTiers = FALSE) {

  if (length(gr) == 0) return(invisible())

  n <- length(gr)
  if (is.null(x0)) x0 <- start(gr)
  if (is.null(x1)) x1 <- end(gr)
  if (is.null(labels)) labels <- rep(NA_character_, n)
  if (is.null(orientation)) orientation <- rep("*", n)
  if (length(shape) == 1) shape <- rep(shape, n)


  # Estimate label width for stacking
  label_widths <- strwidth(labels, cex = label_cex)
  x_left <- ifelse(orientation %in% c("-", "-/-", "-/+", "*"), x0, x0 - label_widths)
  x_right <- ifelse(orientation %in% c("-", "-/-", "-/+", "*"), x1 + label_widths, x1)

  # Sort and assign tiers based on annotation+label space
  ord <- order(x_left)
  x0 <- x0[ord]
  x1 <- x1[ord]
  x_left <- x_left[ord]
  x_right <- x_right[ord]
  labels <- labels[ord]
  orientation <- orientation[ord]
  shape <- shape[ord]
  col <- if (length(col) == 1) rep(col, n) else col[ord]

  tiers <- integer(n)
  current_ends <- c()

  for (i in seq_along(x_left)) {
    placed <- FALSE
    for (t in seq_along(current_ends)) {
      if (x_left[i] > current_ends[t]) {
        current_ends[t] <- x_right[i]
        tiers[i] <- t
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      current_ends <- c(current_ends, x_right[i])
      tiers[i] <- length(current_ends)
    }
  }

  # Use y0s/y1s if provided; otherwise compute from ybase + tiers
  if (is.null(y0s) || is.null(y1s)) {
    y0s <- ybase + (tiers - 1) * tier_height
    y1s <- y0s + tier_height * 0.8
  }

  if (!returnTiers){
    for (i in seq_along(x0)) {
      arrow_width <- 8
      dir <- switch(orientation[i], "+" = 1, "-" = -1, "*" = 1)

      if (shape[i] == "rect") {
        rect(x0[i], y0s[i], x1[i], y1s[i], col = col[i], border = border)
      } else if (shape[i] == "arrow") {
        y_mid <- (y0s[i] + y1s[i]) / 2
        width <- abs(x1[i] - x0[i])

        if (width == 0 || is.na(width)) next  # skip zero-length arrows

        if (orientation[i] == "+" || orientation[i] == "*") {
          arrows(x0[i], y_mid, x1[i], y_mid, length = arrow_width / 200, col = col[i], lwd = 2, angle = 45)
        } else {
          arrows(x1[i], y_mid, x0[i], y_mid, length = arrow_width / 200, col = col[i], lwd = 2, angle = 45)
        }

      } else if (shape[i] == "multi-arrow") {
        y_mid <- (y0s[i] + y1s[i]) / 2
        width <- abs(x1[i] - x0[i])

        if (width == 0 || is.na(width)) next

        # Settings
        arrow_spacing <- 10
        arrowhead_len <- 10

        # Orientation-aware setup
        if (orientation[i] == "+" || orientation[i] == "*") {
          x_start <- x0[i]
          x_end <- x1[i]
          by <- arrow_spacing
          head_from <- function(x) arrows(x - arrowhead_len, y_mid, x, y_mid, length = 0.05, col = col[i], lwd = 1, angle = 45)
        } else {
          x_start <- x1[i]
          x_end <- x0[i]
          by <- -arrow_spacing
          head_from <- function(x) arrows(x + arrowhead_len, y_mid, x, y_mid, length = 0.05, col = col[i], lwd = 1, angle = 45)
        }

        # Draw base line
        segments(x0[i], y_mid, x1[i], y_mid, col = col[i], lwd = 1.5)

        # Internal arrows
        if ((by > 0 && x_end - x_start > 2 * arrow_spacing) ||
            (by < 0 && x_start - x_end > 2 * abs(arrow_spacing))) {
          x_seq <- seq(x_start + by, x_end - by, by = by)
          for (x in x_seq) head_from(x)
        }

        # Terminal arrow
        head_from(x_end)

      } else if (shape[i] == "marker") {
        mid_y <- (y0s[i] + y1s[i]) / 2
        width <- abs(x1[i] - x0[i])

        if (orientation[i] == "+" || orientation[i] == "*") {
          if (width < arrow_width) {
            polygon(
              x = c(x0[i], x1[i], x0[i]),
              y = c(y0s[i], mid_y, y1s[i]),
              col = col[i], border = border
            )
          } else {
            polygon(
              x = c(x0[i], x1[i] - arrow_width, x1[i], x1[i] - arrow_width, x0[i]),
              y = c(y0s[i], y0s[i], mid_y, y1s[i], y1s[i]),
              col = col[i], border = border
            )
          }
        } else {
          if (width < arrow_width) {
            polygon(
              x = c(x1[i], x0[i], x1[i]),
              y = c(y0s[i], mid_y, y1s[i]),
              col = col[i], border = border
            )
          } else {
            polygon(
              x = c(x1[i], x0[i] + arrow_width, x0[i], x0[i] + arrow_width, x1[i]),
              y = c(y0s[i], y0s[i], mid_y, y1s[i], y1s[i]),
              col = col[i], border = border
            )
          }
        }
      }

      # Draw label
      label_side <- if (orientation[i] == "-") "right" else "left"
      label_x <- if (label_side == "left") x0[i] - label_offset else x1[i] + label_offset
      label_adj <- if (label_side == "left") 1 else 0

      if (!is.na(labels[i]) && nzchar(labels[i])) {
        text(x = label_x,
             y = (y0s[i] + y1s[i]) / 2,
             labels = labels[i],
             adj = c(label_adj, 0.5),
             cex = label_cex)
      }
    }
}

  if (returnTiers) {
    return(tiers)
  }
}




## SeqIdeogram ----

#' stainToColor
#'
#' Converts UCSC-style cytoband stains into plot colors.
#' @param stains A character vector of gieStain values.
#' @return A character vector of color hex codes.
#' @export
stainToColor <- function(stains) {
  stain_colors <- c(
    gneg = "whitesmoke",
    gpos25 = "grey75",
    gpos50 = "grey50",
    gpos75 = "grey25",
    gpos100 = "grey0",
    gvar = "#67a1ca",
    stalk = "#f063a7",     # lavender for stalk
    acen = "#f34e4e"       # red for centromeres
  )

  result <- stain_colors[as.character(stains)]
  result[is.na(result)] <- "#CCCCCC"  # fallback gray
  return(result)
}



#' SeqIdeogram Class
#' @description S4 class for ideogram annotations.
#' @export
setClass("SeqIdeogram",
         contains = "SeqAnnotation",
         slots = list(
           gr = "GRanges",
           color = "character",
           y0 = "numeric",
           y1 = "numeric",
           x0 = "numeric",
           x1 = "numeric"
         ))



#' SeqIdeogram Constructor
#'
#' @description Constructs a SeqIdeogram object for drawing cytobands.
#' @param cytobands GRanges object with cytoband data (chrom, start, end, gieStain).
#' @param y0 Numeric, lower y position for ideogram (default 0).
#' @param y1 Numeric, upper y position for ideogram (default 1).
#' @return A SeqIdeogram object.
#' @export
SeqIdeogram <- function(cytobands = NULL, y0 = 0, y1 = 1) {
  if (is.null(cytobands)) {
    cytobands <- loadCytobands()
  }

  # Ensure it's a GRanges object
  if (!is(cytobands, "GRanges")) {
    cytobands <- GRanges(
      seqnames = cytobands$chrom,
      ranges = IRanges(start = cytobands$chromStart, end = cytobands$chromEnd),
      gieStain = cytobands$gieStain
    )
  }

  # Map stains to colors
  stains <- mcols(cytobands)$gieStain
  colors <- stainToColor(stains)

  # Precompute x0/x1 as genomic coordinates
  x0 <- start(cytobands)
  x1 <- end(cytobands)

  # Create SeqIdeogram object
  new("SeqIdeogram",
      gr = cytobands,
      color = colors,
      x0 = as.numeric(x0),
      x1 = as.numeric(x1),
      y0 = rep(y0, length(cytobands)),
      y1 = rep(y1, length(cytobands))
      )
}



#' drawIdeogram
#'
#' @description Draws an ideogram track using cytoband rectangles.
#' @param feat A SeqIdeogram object.
#' @param layout Layout object from setupGlobalLayout().
#' @param track_idx Integer index of the current track.
#' @return None. Draws directly to the active plot.
#' @export
drawIdeogram <- function(feat, layout, track_idx) {
  if (!inherits(feat, "SeqIdeogram")) stop("feat must be a SeqIdeogram")

  gr <- feat@gr
  if (length(gr) == 0) return()

  ov <- findOverlaps(gr, layout$global_windows, select = "first")
  valid <- !is.na(ov)
  if (!any(valid)) return()

  x0_abs <- globalTransformX(feat@x0[valid], ov[valid], layout)
  x1_abs <- globalTransformX(feat@x1[valid], ov[valid], layout)
  t_vec <- rep(track_idx, sum(valid))

  y0_abs <- globalTransformY(feat@y0[valid], t_vec, layout,
                             track_y_min = rep(0, length(t_vec)),
                             track_y_max = rep(1, length(t_vec)))
  y1_abs <- globalTransformY(feat@y1[valid], t_vec, layout,
                             track_y_min = rep(0, length(t_vec)),
                             track_y_max = rep(1, length(t_vec)))

  stain <- mcols(gr)$gieStain[valid]
  color <- feat@color[valid]

  for (i in seq_along(x0_abs)) {
    if (stain[i] == "acen") {
      # Draw centromere triangle (hourglass appearance)
      left_triangle <- FALSE
      right_triangle <- FALSE

      if (i < length(x0_abs) && stain[i + 1] == "acen") {
        left_triangle <- TRUE
      } else if (i > 1 && stain[i - 1] == "acen") {
        right_triangle <- TRUE
      }

      if (left_triangle) {
        polygon(
          x = c(x0_abs[i], x1_abs[i], x0_abs[i]),
          y = c(y0_abs[i], (y0_abs[i] + y1_abs[i]) / 2, y1_abs[i]),
          col = color[i], border = NA, lwd = 0.5
        )
      } else if (right_triangle) {
        polygon(
          x = c(x1_abs[i], x0_abs[i], x1_abs[i]),
          y = c(y0_abs[i], (y0_abs[i] + y1_abs[i]) / 2, y1_abs[i]),
          col = color[i], border = NA, lwd = 0.5
        )
      } else {
        # fallback for isolated centromere
        rect(x0_abs[i], y0_abs[i], x1_abs[i], y1_abs[i], col = color[i], border = NA, lwd = 0.5)
      }

    } else {
      # Standard band
      rect(x0_abs[i], y0_abs[i], x1_abs[i], y1_abs[i], col = color[i], border = NA, lwd = 0.5)
    }
  }
}




#' SeqElement-Class
#' @name SeqElement-Class
#' @description
#' S4 parent class for SeqFeature, SeqLink, and SeqAnnotation
#' @author Andrew Lynch
#' @exportClass SeqElement
setClassUnion("SeqElement", c("SeqFeature", "SeqLink", "SeqAnnotation"))



#' SeqTrack-Class
#' @name SeqTrack-Class
#' @description
#' S4 class for SeqTracks
#' @author Andrew Lynch
#' @exportClass SeqTrack
setClass("SeqTrack",
  slots = list(
    features = "list",      # A list of SeqFeature objects
    seqWindows = "GRanges", # Derived genomic windows for each seqname (from input data)
    secondaryAxis = "logical",  # Whether a second y-axis is used
    xAxis = "logical",      # Whether to draw an x-axis (default FALSE)
    yAxis = "logical",      # Whether to draw a y-axis (default TRUE)
    yTitle = "character",    # Y-axis title
    trackBackgroundColor = "character",
    trackBorderColor = "character",
    windowBackgroundColor = "character",
    windowBorderColor = "character",
    fixedStubs = "logical"
  )
)



## SeqGene ----
#' AnnoGene-Class
#' @description Annotation object for grouped exon visualization
#' @export
setClass("AnnoGene",
         contains = "AnnoFeature",
         slots = list(
           gene_id = "character",
           genes = "character"
         )
)


#' AnnoGene
#' @description Draws grouped gene annotations as exon boxes connected by multi-arrow lines
#' @export
AnnoGene <- function(gr,
                     geneCol = "gene_name",
                     genes = NULL,
                     strandCol = NULL,
                     color = "gray30",
                     colorCol = NULL,
                     shape = "rect"
                     ) {

  if (!inherits(gr, "GRanges")) stop("Input must be a GRanges object.")

  # Filter to exon rows (if present)
  if ("type" %in% names(mcols(gr))) {
    exon_mask <- mcols(gr)$type == "exon"
    if (any(exon_mask)) gr <- gr[exon_mask]
  }

  # Check geneCol exists
  if (!geneCol %in% names(mcols(gr))) {
    stop(paste0("Column '", geneCol, "' not found in metadata."))
  }

  # Apply gene filter if provided
  if (!is.null(genes)) {
    gene_ids_all <- as.character(mcols(gr)[[geneCol]])
    gr <- gr[gene_ids_all %in% genes]
  }

  # Collapse to unique exon ranges per gene
  gene_ids <- as.character(mcols(gr)[[geneCol]])
  gr$.__gene_id__ <- gene_ids
  gr <- gr[!duplicated(paste0(seqnames(gr), ":", start(gr), "-", end(gr), "_", gene_ids))]

  n <- length(gr)
  y0 <- rep(0, n)
  y1 <- rep(1, n)

  # Strand handling
  orientation <- if (!is.null(strandCol)) {
    as.character(mcols(gr)[[strandCol]])
  } else {
    as.character(strand(gr))
  }
  orientation[!orientation %in% c("+", "-")] <- "*"

  # Default labels = gene name
  labels <- gene_ids
  labels[is.na(labels)] <- ""

  # Colors
  colors <- if (!is.null(colorCol) && colorCol %in% names(mcols(gr))) {
    as.character(mcols(gr)[[colorCol]])
  } else {
    rep(color, n)
  }

  # Shape
  shapes <- rep(shape, n)

  new("AnnoGene",
      gr = gr,
      labels = labels,
      color = colors,
      shape = shapes,
      orientation = orientation,
      gene_id = gene_ids,
      genes = unique(gene_ids),
      y0 = y0,
      y1 = y1)
}



#' drawAnnoGene
#' @description Draws grouped gene annotations as exon boxes connected by multi-arrow lines
#' @export
#' drawAnnoGene
#' @description Draws grouped gene annotations as exon boxes connected by multi-arrow lines
#' @export
drawAnnoGene <- function(ag, layout, track_idx,
                         track_y_min = 0, track_y_max = 1,
                         totalTracks,
                         arrow_spacing = 10,
                         arrowhead_len = 1,
                         exon_height_ratio = 0.8,
                         label_cex = 0.7,
                         label_offset = 5) {

  if (!inherits(ag, "AnnoGene")) stop("Input must be an AnnoGene object")
  gr <- ag@gr
  if (length(gr) == 0) return(invisible())

  # Pre-filter to exons that actually overlap visible windows
  ov <- findOverlaps(gr, layout$global_windows, select = "first")
  valid <- !is.na(ov)
  if (!any(valid)) return(invisible())

  # Keep only overlapping exons and metadata
  gr           <- gr[valid]
  gene_ids     <- ag@gene_id[valid]
  labels       <- ag@labels[valid]
  colors       <- ag@color[valid]
  orientations <- ag@orientation[valid]

  # Group by gene
  gene_map        <- split(gr, gene_ids)
  label_map       <- tapply(labels, gene_ids, function(x) x[1])
  color_map       <- tapply(colors, gene_ids, function(x) x[1])
  orientation_map <- tapply(orientations, gene_ids, function(x) x[1])

  label_pad <- 50000

  gene_ranges <- GRanges(
    seqnames = sapply(gene_map, function(g) as.character(seqnames(g)[1])),
    ranges = IRanges(
      start = mapply(function(g, ori) {
        if (ori == "+") min(start(g)) - label_pad else min(start(g))
      }, gene_map, orientation_map[names(gene_map)]),
      end = mapply(function(g, ori) {
        if (ori == "-") max(end(g)) + label_pad else max(end(g))
      }, gene_map, orientation_map[names(gene_map)])
    ),
    gene_id = names(gene_map)
  )


  assignBrickTiers <- function(gr) {
    gr <- gr[order(seqnames(gr), start(gr))]
    tiers <- integer(length(gr))
    ranges_by_seq <- split(gr, seqnames(gr))

    for (seq in names(ranges_by_seq)) {
      group <- ranges_by_seq[[seq]]
      end_last <- c()
      for (i in seq_along(group)) {
        s <- start(group)[i]
        tier <- which(s >= end_last)
        if (length(tier) == 0) {
          tiers[which(gr == group[i])] <- length(end_last) + 1
          end_last <- c(end_last, end(group)[i])
        } else {
          t <- tier[1]
          tiers[which(gr == group[i])] <- t
          end_last[t] <- end(group)[i]
        }
      }
    }

    return(tiers)
  }

  assignBrickTiersAbsolute <- function(x0, x1) {
    o <- order(x1)
    x0 <- x0[o]
    x1 <- x1[o]
    n <- length(x0)
    tiers <- integer(n)
    end_last <- numeric(0)

    for (i in seq_len(n)) {
      s <- x0[i]
      tier <- which(s >= end_last)
      if (length(tier) == 0) {
        t <- length(end_last) + 1
        end_last <- c(end_last, x1[i])
      } else {
        t <- tier[1]
        end_last[t] <- x1[i]
      }
      tiers[o[i]] <- t
    }

    return(tiers)
  }

  # Compute one gene span per group
  gene_ids <- names(gene_map)

  gene_info <- lapply(gene_ids, function(gid) {
    g <- gene_map[[gid]]
    list(
      gid = gid,
      seqname = as.character(seqnames(g)[1]),
      start = min(start(g)),
      end = max(end(g)),
      strand = orientation_map[[gid]]
    )
  })
  gene_info <- do.call(rbind, lapply(gene_info, as.data.frame))

  # Transform to absolute x-coordinates
  gene_ov <- findOverlaps(GRanges(gene_info$seqname, IRanges(gene_info$start, gene_info$end)),
                          layout$global_windows, select = "first")

  x0_abs <- globalTransformX(gene_info$start, gene_ov, layout)
  x1_abs <- globalTransformX(gene_info$end, gene_ov, layout)

  # Estimate label padding in absolute space
  label_pad <- (layout$x1 - layout$x0) * 0.03  # ~3% of width

  # Adjust for strand
  x0_label <- ifelse(gene_info$strand == "+", x0_abs - label_pad, x0_abs)
  x1_label <- ifelse(gene_info$strand == "-", x1_abs + label_pad, x1_abs)

  # Assign non-overlapping tiers
  gene_tiers <- assignBrickTiersAbsolute(x0_label, x1_label)
  names(gene_tiers) <- gene_ids
  n_tiers <- max(gene_tiers)

  tier_height <- (track_y_max - track_y_min) / n_tiers

  for (i in seq_along(gene_ids)) {
    gid <- names(gene_map)[i]
    gene_exons <- gene_map[[i]]
    if (length(gene_exons) == 0) next

    label <- if (!is.null(label_map[[gid]])) label_map[[gid]] else gid
    color <- if (!is.null(color_map[[gid]])) color_map[[gid]] else "gray30"
    orientation <- if (!is.null(orientation_map[[gid]])) orientation_map[[gid]] else "*"

    # Y-coordinates (one line per gene)
    tier_idx <- gene_tiers[[gid]]
    exon_y0 <- track_y_min + (tier_idx - 1) * tier_height
    exon_y1 <- exon_y0 + tier_height * exon_height_ratio
    y_mid   <- (exon_y0 + exon_y1) / 2

    # Map exon coords to plot space
    exon_ov <- findOverlaps(gene_exons, layout$global_windows, select = "first")
    exon_valid <- !is.na(exon_ov)
    if (!any(exon_valid)) next

    gene_exons <- gene_exons[exon_valid]
    exon_ov <- exon_ov[exon_valid]
    x0_abs <- globalTransformX(start(gene_exons), exon_ov, layout)
    x1_abs <- globalTransformX(end(gene_exons), exon_ov, layout)

    t_scalar <- track_idx

    # Transform Y coordinates into absolute plotting space
    exon_y0_abs <- globalTransformY(exon_y0, t_scalar, layout,
                                    track_y_min = rep(track_y_min, totalTracks),
                                    track_y_max = rep(track_y_max, totalTracks))
    exon_y1_abs <- globalTransformY(exon_y1, t_scalar, layout,
                                    track_y_min = rep(track_y_min, totalTracks),
                                    track_y_max = rep(track_y_max, totalTracks))
    y_mid_abs   <- (exon_y0_abs + exon_y1_abs) / 2

    # Gene span (backbone)
    gene_start <- min(start(gene_exons))
    gene_end   <- max(end(gene_exons))
    gene_ov    <- exon_ov[which.min(start(gene_exons))]
    if (is.na(gene_ov)) next

    line_x0 <- globalTransformX(gene_start, gene_ov, layout)
    line_x1 <- globalTransformX(gene_end, gene_ov, layout)

    #print(c(line_x0 = line_x0, line_x1 = line_x1, exon_y0_abs = exon_y0_abs, exon_y1_abs = exon_y1_abs)) #!

    # Orientation-aware arrow logic
    if (orientation == "+" || orientation == "*") {
      x_start <- line_x0 + arrow_spacing
      x_end   <- line_x1 - arrow_spacing
      by      <- arrow_spacing
      draw_terminal_arrow <- function() {
        arrows(line_x1 - arrowhead_len, y_mid_abs, line_x1, y_mid_abs, length = 0.05, col = color, lwd = 1, angle = 45)
      }
    } else {
      x_start <- line_x1 - arrow_spacing
      x_end   <- line_x0 + arrow_spacing
      by      <- -arrow_spacing
      draw_terminal_arrow <- function() {
        arrows(line_x0 + arrowhead_len, y_mid_abs, line_x0, y_mid_abs, length = 0.05, col = color, lwd = 1, angle = 45)
      }
    }

    # Draw gene backbone
    segments(line_x0, y_mid_abs, line_x1, y_mid_abs, col = color, lwd = 1.5)

    # if ((by > 0 && x_start < x_end) || (by < 0 && x_start > x_end)) {
    #   x_seq <- seq(x_start, x_end, by = by)
    #   for (x in x_seq) {
    #     if (orientation == "+" || orientation == "*") {
    #       arrows(x - arrowhead_len, y_mid_abs, x, y_mid_abs, length = 0.05, col = color, lwd = 1, angle = 45)
    #     } else {
    #       arrows(x + arrowhead_len, y_mid_abs, x, y_mid_abs, length = 0.05, col = color, lwd = 1, angle = 45)
    #     }
    #   }
    # }

    draw_terminal_arrow()

    # Draw exon rectangles
    for (j in seq_along(x0_abs)) {
      rect(x0_abs[j], exon_y0_abs, x1_abs[j], exon_y1_abs, col = color, border = color)
    }

    # Draw label
    if (!is.na(label) && nzchar(label)) {
      label_x   <- if (orientation == "+") line_x0 - label_offset else line_x1 + label_offset
      label_adj <- if (orientation == "+") 1 else 0
      text(x = label_x, y = y_mid_abs, labels = label, adj = c(label_adj, 0.5),
           cex = label_cex)
    }
  }

  invisible()
}






# SeqTrack ----
#' SeqTrack
#' @name SeqTrack
#' @description
#' Combines multiple SeqFeature, SeqLink, or SeqAnnotation objects into a track.
#'
#' @param featureList List of feature objects.
#' @param secondaryAxis Logical; draw a second y-axis.
#' @param xAxis Logical; draw an x-axis.
#' @param yAxis Logical; draw a y-axis.
#' @param yTitle Character; title of the y-axis.
#' @return A SeqTrack object.
#' @export
SeqTrack <- function(featureList, secondaryAxis = FALSE, xAxis = FALSE, yAxis = TRUE, yTitle = "", trackBackgroundColor = "transparent", trackBorderColor = "transparent", windowBackgroundColor = "transparent", windowBorderColor = "grey50", fixedStubs = T) {
  # valid_classes <- c("SeqFeature", "SeqPoint", "SeqSegment", "SeqBar", "SeqLink", "SeqAnnotation", "SeqIdeogram")
  #
  # if (!all(sapply(featureList, function(x) class(x)[1] %in% valid_classes))) {
  #   stop("All elements in featureList must be one of: ", paste(valid_classes, collapse = ", "))
  # }
  if (!all(sapply(featureList, function(x) {
    is(x, "SeqFeature") || is(x, "SeqLink") || is(x, "SeqAnnotation")
  }))) {
    stop("All elements in featureList must be SeqFeature, SeqLink, or SeqAnnotation objects.")
  }

  grAll <- do.call(c, lapply(featureList, function(sf) {
    if (is(sf, "SeqFeature")) {
      return(sf@gr)
    } else if (is(sf, "SeqLink")) {
      gr1 <- sf@gr1
      gr2 <- sf@gr2
      if (!is(gr1, "GRanges") || !is(gr2, "GRanges")) {
        stop("gr1 or gr2 is not a GRanges object in a SeqLink")
      }
      return(c(gr1, gr2))
    } else if (is(sf, "SeqAnnotation")) {
      return(GenomicRanges::GRanges())  # annotations don’t define plotting extent
    } else {
      stop("Unknown feature type in featureList")
    }
  }))

  grWindows <- sort(range(GenomicRanges::GRangesList(split(grAll, seqnames(grAll)))))

  new("SeqTrack",
      features = featureList,
      seqWindows = unlist(grWindows),
      secondaryAxis = secondaryAxis,
      xAxis = xAxis,
      yAxis = yAxis,
      yTitle = yTitle,
      trackBackgroundColor = trackBackgroundColor,
      trackBorderColor = trackBorderColor,
      windowBackgroundColor = windowBackgroundColor,
      windowBorderColor = windowBorderColor,
      fixedStubs = fixedStubs)
}


# SeqPlot ----
## SeqPlot Class ----

#' SeqPlot-Class
#'
#' @name SeqPlot-Class
#' @description
#' S4 class for SeqPlots
#' @author Andrew Lynch
#' @exportClass SeqPlot
setClass("SeqPlot",
  slots = list(
    tracks           = "list",     # A list of SeqTrack objects
    hideEmptyWindows = "logical",  # Option to hide unused windows in each track
    trackGap         = "numeric",  # Gap between tracks
    windowGap        = "numeric",  # Gap between genome windows (x-axis)
    margin           = "list",     # Margin for the entire plot (list with top, right, bottom, left)
    trackHeights     = "numeric"   # Relative heights for each track
  )
)

## SeqPlot Constructor ----
#' Create a SeqPlot object
#'
#' Stacks multiple SeqTracks and manages layout, margins, and window order.
#'
#' @param trackList List of SeqTrack objects.
#' @param hideEmptyWindows Logical; hide windows with no data.
#' @param trackGap Numeric; vertical gap between tracks.
#' @param windowGap Numeric; horizontal gap between genome windows.
#' @param margin List with plot margins (top, right, bottom, left).
#' @param trackHeights Optional vector of relative heights for each track.
#' @return A SeqPlot object.
#' @export
SeqPlot <- function(trackList, hideEmptyWindows = FALSE, trackGap = 10, windowGap = 20,
                    margin = list(top = 20, right = 10, bottom = 50, left = 100),
                    trackHeights = NULL) {

  if (!all(sapply(trackList, function(x) is(x, "SeqTrack")))) {
    stop("All elements in trackList must be SeqTrack objects.")
  }

  totalTracks <- length(trackList)

  if (is.null(trackHeights)) {
    trackHeights <- rep(1, totalTracks)  # default: equal height
  }

  if (length(trackHeights) != totalTracks) {
    stop("Length of trackHeights must match number of tracks.")
  }

  # Normalize heights so they represent proportions
  trackHeights <- trackHeights / sum(trackHeights)

  new("SeqPlot",
      tracks = trackList,
      hideEmptyWindows = hideEmptyWindows,
      trackGap = trackGap,
      windowGap = windowGap,
      margin = margin,
      trackHeights = trackHeights)
}

## SeqPlot Method ----

setGeneric("plotSeqPlot", function(sp, globalWindows = NULL, windowOrder = NULL)
  standardGeneric("plotSeqPlot"))

setMethod("plotSeqPlot", "SeqPlot", function(sp, globalWindows = NULL, windowOrder = NULL) {

  # Use globalWindows from argument, sp, or defaultGenomeWindows()
  if (is.null(globalWindows)) {
    global_windows <- defaultGenomeWindows()
  } else {
    global_windows <- globalWindows
  }

  # If windowOrder is provided, reorder global_windows accordingly
  if (!is.null(windowOrder)) {
    global_windows <- global_windows[as.character(seqnames(global_windows)) %in% windowOrder]
    global_windows <- global_windows[order(factor(as.character(seqnames(global_windows)), levels = windowOrder), start(global_windows))]
  }

  ## Set up the layout using the new module (returns a list with global_* parameters)
  layout_info <- setupGlobalLayout(sp, globalWindows = globalWindows, plotWidth = 1000)
  margin <- layout_info$global_margin
  global_window_starts <- layout_info$global_window_starts
  global_window_ends <- layout_info$global_window_ends
  global_track_bottoms <- layout_info$global_track_bottoms
  global_track_tops <- layout_info$global_track_tops
  totalWidth <- layout_info$global_plot_width
  totalHeight <- layout_info$global_plot_height  # includes margins

  # Define the plotting origin
  global_origin_x <- margin$left
  global_origin_y <- margin$bottom

  # Initialize plot canvas
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, totalWidth), ylim = c(0, totalHeight),
       type = "n", axes = FALSE, xlab = "", ylab = "")

  # Draw background for each track
  totalTracks <- length(sp@tracks)
  for (i in seq_len(totalTracks)) {
    track <- sp@tracks[[i]]
    rect(global_origin_x, global_track_bottoms[i],
         global_origin_x + 1000, global_track_tops[i],
         col = track@trackBackgroundColor, border = track@trackBorderColor)
  }

  # --- Draw any SeqIdeogram objects ---
  for (track_idx in seq_len(totalTracks)) {
    track <- sp@tracks[[track_idx]]
    for (sf in track@features) {
      if (inherits(sf, "SeqIdeogram")) {
        drawIdeogram(sf, layout_info, track_idx)
      }
    }
  }

  # Optional: Draw background box for each window inside each track
  for (track_idx in seq_len(totalTracks)) {
    track <- sp@tracks[[track_idx]]
    for (win_idx in seq_along(layout_info$global_windows)) {
      rect(layout_info$global_window_starts[win_idx],
           layout_info$global_track_bottoms[track_idx],
           layout_info$global_window_ends[win_idx],
           layout_info$global_track_tops[track_idx],
           col = track@windowBackgroundColor, border = track@windowBorderColor)
    }
  }

  ## Process each track for standard (non-spanning) features:
  for (track_idx in seq_len(totalTracks)) {
    track <- sp@tracks[[track_idx]]
    # Compute track-level y range from features
    track_y_vals <- unlist(lapply(track@features, function(sf) {
      yvals <- NULL

      # Use y0/y1 if present (for all SeqElements)
      if (is(sf, "SeqElement")) {
        yvals <- c(yvals, sf@y0, sf@y1)
      }

      # Include height only for SeqArch and SeqRecon
      if (inherits(sf, "SeqArch") || inherits(sf, "SeqRecon")) {
        yvals <- c(yvals, sf@height)
      }

      # Handle NA and infinite cleanup
      yvals[is.infinite(yvals)] <- NA

      # Fallback for annotations like AnnoGene
      if (inherits(sf, "AnnoGene")) {
        yvals <- c(yvals, 0, 1)
      }

      return(yvals)
    }))

    if (length(track_y_vals) == 0) {
      track_y_min <- 0
      track_y_max <- 1
    } else {
      track_y_min <- min(track_y_vals, na.rm = TRUE)
      track_y_max <- max(track_y_vals, na.rm = TRUE)
    }
    if (!is.finite(track_y_min) || !is.finite(track_y_max) || track_y_min == track_y_max) {
      track_y_min <- 0
      track_y_max <- 1
    }

    this_track_bottom <- global_track_bottoms[track_idx]
    this_track_top <- global_track_tops[track_idx]
    thisTrackHeight <- this_track_top - this_track_bottom

    # Containers for vectorized drawing of standard features
    seg_x0 <- numeric(); seg_y0 <- numeric(); seg_x1 <- numeric(); seg_y1 <- numeric(); seg_col <- character()
    pt_x <- numeric(); pt_y <- numeric(); pt_col <- character()
    box_x0 <- numeric(); box_y0 <- numeric(); box_x1 <- numeric(); box_y1 <- numeric(); box_col <- character()

    # Loop over annotation features
    for (sf in track@features) {
      if (inherits(sf, "AnnoGene")) {
        drawAnnoGene(sf,
                     layout = layout_info,
                     track_idx = track_idx,
                     track_y_min = track_y_min,
                     track_y_max = track_y_max,
                     totalTracks = length(sp@tracks))

      } else if (is(sf, "AnnoFeature")) {

        gr <- sf@gr
        ov <- findOverlaps(gr, layout_info$global_windows, select = "first")
        valid <- !is.na(ov)
        if (!any(valid)) next

        gr <- gr[valid]
        win_idx <- ov[valid]
        labels <- sf@labels[valid]
        strand_vals <- as.character(strand(gr))
        orientation <- ifelse(strand_vals %in% c("+", "-"), strand_vals, "*")

        shape <- if (!is.null(sf@shape)) sf@shape[valid] else rep("rect", length(gr))
        col <- sf@color[valid]

        x0_abs <- globalTransformX(start(gr), win_idx, layout_info)
        x1_abs <- globalTransformX(end(gr), win_idx, layout_info)

        # For each annotation, assign it to this track index
        t_vec <- rep(track_idx, length(gr))

        # Compute y0/y1 for each tier
        # Compute y0/y1 for each tier (before transform)
        tiers <- drawAnnoFeature(
          gr = gr,
          labels = labels,
          orientation = orientation,
          shape = shape,
          x0 = x0_abs,
          x1 = x1_abs,
          col = col,
          returnTiers = TRUE
        )

        n_tiers <- max(tiers)
        y0_rel <- track_y_min + (tiers - 1) * (track_y_max - track_y_min) / n_tiers
        y1_rel <- y0_rel + (track_y_max - track_y_min) / n_tiers * 0.8

        # Convert to absolute plotting coordinates
        y0_abs <- globalTransformY(y0_rel, t_vec, layout_info,
                                   track_y_min = track_y_min[track_idx],
                                   track_y_max = track_y_max[track_idx])
        y1_abs <- globalTransformY(y1_rel, t_vec, layout_info,
                                   track_y_min = track_y_min[track_idx],
                                   track_y_max = track_y_max[track_idx])

        # Optionally: overwrite x0/x1 in GRanges to preserve for labels, if needed
        # mcols(gr)$x0_abs <- x0_abs
        # mcols(gr)$x1_abs <- x1_abs

        # Replace the GRanges' ranges with transformed values (tricky but works)
        ranges(gr) <- IRanges(start = x0_abs, end = x1_abs)
        drawAnnoFeature(
          gr = gr,
          labels = labels,
          orientation = orientation,
          shape = shape,
          x0 = x0_abs,
          x1 = x1_abs,
          y0s = y0_abs,
          y1s = y1_abs,
          col = col
        )

      }
    }

    # Loop over each feature that is not a spanning (SeqLink) feature.
    for (sf in track@features) {
      if (inherits(sf, "SeqLink")) next

      gr <- sf@gr
      ov <- findOverlaps(gr, global_windows, select = "first")
      valid <- !is.na(ov)

      if (!any(valid)) next
      win_idx_vec <- ov[valid]

      # Transform x coordinates using our globalTransformX helper
      x0_abs <- globalTransformX(sf@x0[valid], win_idx_vec, layout_info)
      x1_abs <- globalTransformX(sf@x1[valid], win_idx_vec, layout_info)
      # Transform y coordinates using globalTransformY
      # Note: For non-link features, the track index is simply track_idx.
      t_vec <- rep(track_idx, sum(valid))
      y0_abs <- globalTransformY(sf@y0[valid], t_vec, layout_info,
                                 track_y_min = rep(track_y_min, totalTracks),
                                 track_y_max = rep(track_y_max, totalTracks))
      y1_abs <- globalTransformY(sf@y1[valid], t_vec, layout_info,
                                 track_y_min = rep(track_y_min, totalTracks),
                                 track_y_max = rep(track_y_max, totalTracks))

      if (inherits(sf, "SeqSegment")) {
        seg_x0 <- c(seg_x0, x0_abs)
        seg_y0 <- c(seg_y0, y0_abs)
        seg_x1 <- c(seg_x1, x1_abs)
        seg_y1 <- c(seg_y1, y1_abs)
        seg_col <- c(seg_col, sf@color[valid])
      } else if (inherits(sf, "SeqPoint")) {
        pt_x <- c(pt_x, x0_abs)
        pt_y <- c(pt_y, y0_abs)
        pt_col <- c(pt_col, sf@color[valid])
      } else if (inherits(sf, "SeqBar")) {
        box_x0 <- c(box_x0, x0_abs)
        box_y0 <- c(box_y0, y0_abs)
        box_x1 <- c(box_x1, x1_abs)
        box_y1 <- c(box_y1, y1_abs)
        box_col <- c(box_col, sf@color[valid])
      }

    }

    # Draw vectorized non-link features
    if (length(seg_x0) > 0) segments(seg_x0, seg_y0, seg_x1, seg_y1, col = seg_col, lwd = 2)
    if (length(pt_x) > 0) points(pt_x, pt_y, pch = 16, col = pt_col, cex = 0.5)
    if (length(box_x0) > 0) rect(box_x0, box_y0, box_x1, box_y1, col = box_col, border = NA)

    # Draw x-axis below each window if enabled for the track
    if (track@xAxis) {
      for (win_idx in seq_len(length(layout_info$global_windows))) {
        x0 <- layout_info$global_window_starts[win_idx]
        x1 <- layout_info$global_window_ends[win_idx]
        chrom <- gsub("^chr", "", as.character(seqnames(layout_info$global_windows)[win_idx]))
        scaled_start <- start(layout_info$global_windows)[win_idx] * layout_info$global_scale_factors[win_idx]
        scaled_end <- end(layout_info$global_windows)[win_idx] * layout_info$global_scale_factors[win_idx]

        # Determine unit label based on scale factor
        sf <- layout_info$global_scale_factors[win_idx]
        if (abs(sf - 1e-6) < 1e-9) {
          unit_label <- "Mb"
        } else if (abs(sf - 1e-3) < 1e-9) {
          unit_label <- "kb"
        } else if (abs(sf - 1) < 1e-9) {
          unit_label <- "bp"
        } else {
          exponent <- round(log10(sf), 1)
          unit_label <- paste0("e", exponent, " bp")
        }

        tick_length <- 4
        label_cex <- 0.75
        axis_col <- "black"

        # X-axis line and ticks
        segments(x0, this_track_bottom, x1, this_track_bottom, col = axis_col)
        segments(x0, this_track_bottom, x0, this_track_bottom - tick_length, col = axis_col)
        segments(x1, this_track_bottom, x1, this_track_bottom - tick_length, col = axis_col)

        # Labels
        text(x0, this_track_bottom - tick_length - 2, labels = paste0(round(scaled_start, 1), " ", unit_label),
             cex = label_cex, srt = 90, adj = c(1, 1))
        text(x1, this_track_bottom - tick_length - 2, labels = paste0(round(scaled_end, 1), " ", unit_label),
             cex = label_cex, srt = 90, adj = c(1, 0))
        text((x0 + x1) / 2, this_track_bottom - tick_length - 2, labels = chrom,
             cex = label_cex, pos = 1)
      }
    }

    if (track@yAxis) {
      x_pos <- layout_info$global_window_starts[1]
      tick_length <- 4
      label_cex <- 0.75
      axis_col <- "black"

      y0 <- this_track_bottom
      y1 <- layout_info$global_track_tops[track_idx]

      # Main axis line
      segments(x_pos, y0, x_pos, y1, col = axis_col)

      # Ticks at bottom and top
      segments(x_pos, y0, x_pos - tick_length, y0, col = axis_col)
      segments(x_pos, y1, x_pos - tick_length, y1, col = axis_col)

      # Labels
      text(x_pos - tick_length - 2, y0, labels = round(track_y_min, 2), pos = 2, cex = label_cex)
      text(x_pos - tick_length - 2, y1, labels = round(track_y_max, 2), pos = 2, cex = label_cex)

      # Optional: Y-axis title
      if (nchar(track@yTitle) > 0) {
        title_x <- x_pos - tick_length - 40
        title_y <- (y0 + y1) / 2
        text(title_x, title_y, labels = track@yTitle, font = 1, adj = c(1, 0.5), cex = label_cex)
      }
    }
  }

  ## Process spanning features (SeqLink objects)
  for (track_idx in seq_len(totalTracks)) {
    track <- sp@tracks[[track_idx]]
    for (sf in track@features) {
      if (!inherits(sf, "SeqLink")) next

      # Use findOverlaps to map endpoints to windows
      N <- length(sf@gr1)
      ov1 <- findOverlaps(sf@gr1, global_windows, select = "first")
      ov2 <- findOverlaps(sf@gr2, global_windows, select = "first")

      for (i in seq_len(N)) {
        visible_start <- !is.na(ov1[i])
        visible_end   <- !is.na(ov2[i])

        if (!visible_start && !visible_end) next  # completely off screen — skip

        pt1 <- sf@gr1[i]
        pt2 <- sf@gr2[i]

        x0_abs <- if (visible_start) globalTransformX(start(pt1), ov1[i], layout_info) else NA
        x1_abs <- if (visible_end)   globalTransformX(start(pt2), ov2[i], layout_info) else NA

        track_index_0 <- if (!is.na(sf@t0[i]) && sf@t0[i] != 0) sf@t0[i] else track_idx
        track_index_1 <- if (!is.na(sf@t1[i]) && sf@t1[i] != 0) sf@t1[i] else track_idx

        y0_val <- sf@y0[i]
        y1_val <- sf@y1[i]

        y0_abs <- if (!is.na(y0_val) && visible_start) globalTransformY(y0_val, track_index_0, layout_info,
                                                                        track_y_min = rep(track_y_min, totalTracks),
                                                                        track_y_max = rep(track_y_max, totalTracks)) else NA

        y1_abs <- if (!is.na(y1_val) && visible_end) globalTransformY(y1_val, track_index_1, layout_info,
                                                                      track_y_min = rep(track_y_min, totalTracks),
                                                                      track_y_max = rep(track_y_max, totalTracks)) else NA


        # Top positions for SeqArch types
        if (inherits(sf, "SeqString")) {
          print(paste(y0_abs, y1_abs))
          drawSeqString(
            x0 = x0_abs, y0 = y0_abs,
            x1 = x1_abs, y1 = y1_abs,
            orientation = sf@orientation[i],
            col = sf@color[i],
            lwd = 1.5
          )
        } else if (inherits(sf, "SeqArch") || inherits(sf, "SeqRecon")) {
          height <- sf@height[i]
          top0_abs <- if (visible_start) globalTransformY(height, track_idx, layout_info,
                                                          track_y_min = rep(track_y_min, totalTracks),
                                                          track_y_max = rep(track_y_max, totalTracks)) else NA
          top1_abs <- if (visible_end) globalTransformY(height, track_idx, layout_info,
                                                        track_y_min = rep(track_y_min, totalTracks),
                                                        track_y_max = rep(track_y_max, totalTracks)) else NA
          # Dispatch drawing
          if (inherits(sf, "SeqRecon")) {
            if (visible_start && visible_end) {
              drawSeqArch(
                x0 = x0_abs, y0 = y0_abs,
                x1 = x1_abs, y1 = y1_abs,
                top0 = top0_abs, top1 = top1_abs,
                orientation = sf@orientation[i],
                curve = sf@curve[i],
                col = sf@color[i],
                lwd = 1.2
              )
            } else if (visible_start && !visible_end) {
              track_top <- layout_info$global_track_tops[track_idx]
              drawLinkStub(
                x = x0_abs,
                y = y0_abs,
                direction = "right",
                layout = layout_info,
                track_idx = track_idx,
                height = sf@height[i],
                fixedStubs = track@fixedStubs,
                track_y_min = track_y_min,
                track_y_max = track_y_max,
                col = sf@color[i],
                lwd = 1.2
              )

            } else if (!visible_start && visible_end) {
              track_top <- layout_info$global_track_tops[track_idx]
              drawLinkStub(
                x = x1_abs,
                y = y1_abs,
                direction = "left",
                layout = layout_info,
                track_idx = track_idx,
                height = sf@height[i],
                fixedStubs = track@fixedStubs,
                track_y_min = track_y_min,
                track_y_max = track_y_max,
                col = sf@color[i],
                lwd = 1.2
              )
            }
          } else if (is(sf, "SeqArch") && !is(sf, "SeqRecon")) {
            arch_height <- sf@height[i]
            arch_orientation <- if (!is.na(sf@orientation[i])) sf@orientation[i] else "*"

            track_index_0 <- if (!is.na(sf@t0[i]) && sf@t0[i] != 0) sf@t0[i] else track_idx
            track_index_1 <- if (!is.na(sf@t1[i]) && sf@t1[i] != 0) sf@t1[i] else track_idx

            # Overlap detection for gr1/gr2
            visible_start <- !is.na(ov1[i])
            visible_end   <- !is.na(ov2[i])

            x0_abs <- if (visible_start) globalTransformX(start(sf@gr1[i]), ov1[i], layout_info) else NA
            x1_abs <- if (visible_end)   globalTransformX(start(sf@gr2[i]), ov2[i], layout_info) else NA

            y0_val <- sf@y0[i]
            y1_val <- sf@y1[i]


            y0_abs <- if (visible_start) globalTransformY(y0_val, track_index_0, layout_info,
                                                          track_y_min = rep(0, totalTracks),
                                                          track_y_max = rep(1, totalTracks)) else NA

            y1_abs <- if (visible_end) globalTransformY(y1_val, track_index_1, layout_info,
                                                        track_y_min = rep(0, totalTracks),
                                                        track_y_max = rep(1, totalTracks)) else NA

            top0_abs <- if (visible_start) globalTransformY(arch_height, track_idx, layout_info,
                                                            track_y_min = rep(track_y_min, totalTracks),
                                                            track_y_max = rep(track_y_max, totalTracks)) else NA

            top1_abs <- if (visible_end) globalTransformY(arch_height, track_idx, layout_info,
                                                          track_y_min = rep(track_y_min, totalTracks),
                                                          track_y_max = rep(track_y_max, totalTracks)) else NA

            if (visible_start && visible_end) {
              drawSeqArch(
                x0 = x0_abs, y0 = y0_abs,
                x1 = x1_abs, y1 = y1_abs,
                top0 = top0_abs, top1 = top1_abs,
                orientation = sf@orientation[i],
                curve = sf@curve[i],
                col = sf@color[i],
                lwd = 1.2
              )
            } else if (visible_start && !visible_end) {
              drawLinkStub(
                x = x0_abs,
                y = y0_abs,
                direction = "right",
                layout = layout_info,
                track_idx = track_idx,
                height = sf@height[i],
                fixedStubs = track@fixedStubs,
                track_y_min = track_y_min,
                track_y_max = track_y_max,
                col = sf@color[i],
                lwd = 1.2
              )

            } else if (!visible_start && visible_end) {
              drawLinkStub(
                x = x1_abs,
                y = y1_abs,
                direction = "left",
                layout = layout_info,
                track_idx = track_idx,
                height = sf@height[i],
                fixedStubs = track@fixedStubs,
                track_y_min = track_y_min,
                track_y_max = track_y_max,
                col = sf@color[i],
                lwd = 1.2
              )
            }
          }
        } else if (inherits(sf, "SeqLine")) {
          segments(x0_abs, y0_abs, x1_abs, y1_abs, col = sf@color[i], lwd = 1.2)

        } else if (inherits(sf, "SeqArc")) {
          orientation <- if (!is.na(sf@orientation[i])) sf@orientation[i] else "*"
          drawArc(x0_abs, y0_abs, x1_abs, y1_abs, orientation, col = sf@color[i], lwd = 1.2)

        } else {
          print(class(sf))
          warning("Unknown SeqLink subclass — skipping.")
        }

      }
    }
  }
})


# TODOS ----

# SeqFeatures
# TODO - [X] SeqFeature - Point
# TODO - [X] SeqFeature - Segment
# TODO - [X] SeqFeature - Bar
# TODO - [ ] SeqFeature - Stacked Bar
# TODO - [ ] SeqFeature - Line (path through observations)
# TODO - [ ] SeqFeature - Area (basically line but a filled polygon)
# TODO - [ ] SeqFeature - Stacked Area
# TODO - [ ] SeqFeature - Segment
# TODO - [ ] SeqFeature - Tile (Heatmap with ranges on X and categorical on Y.)
# TODO - [ ] SeqFeature - Text
# TODO - [ ] Plot 'stubs' of data that partially resides outside the plotting window.
# TODO - [ ] Grouping of input observations (most important for line, area, and barplot. Will enable stacking.)

# SeqLinks (see genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1895-9 for nice examples)
# TODO - [X] SeqLink - Arc
# TODO - [X] SeqLink - Arch (optional height argument for y maximum, differs from arch in that it is bounded by vertical lines and added height functionality)
# TODO - [/] SeqLink - ReCon (Special arch that separates types of rearrangements. Based on doi.org/10.1093/bioinformatics/btad719. Places inversions on one track and dup/del on another.
       # - [ ]  Still need to add proper tier lines.
# TODO - [/] SeqLink - String (optional height argument for y maximum)
# TODO - [ ] SeqLink - Barbell (intended to visualize paired reads)
# TODO - [ ] SeqLink - Band (for synteny-like connections)
# TODO - [/] Plot 'stubs' of data whose partner is outside the plotting window.
       # - [ ] Still need to fix stubbed SeqArch height scaling and add stubs to SeqStrings

# SeqAnnotations
# TODO - [X] SeqAnnotation - Ideogram
# TODO - [X] SeqAnnotation - Feature (generic box or arrow marker)
# TODO - [X] SeqAnnotation - Genes/Exons/UTR/ETC...
# TODO - [ ] SeqAnnotation - Transcripts
# TODO - [ ] SeqAnnotation - Enhancers
# TODO - [ ] SeqAnnotation - Legend
# TODO - [ ] SeqAnnotation - Alignment
# TODO - [ ] SeqAnnotation - Zoom
# TODO - [ ] SeqAnnotation - Text

# Window/Track Control
# TODO - [ ] Manually set Y axes per track
# TODO - [ ] Convert to grid graphics (rather than single canvas)
# TODO - [ ] Add functionality for customWindows at the Seq Track level. I.e., not all tracks have the same windows.
# TODO - [ ] Individual track and window gap control
# TODO - [ ] Fix Y axis scaling. Track currently scales to maximum Y even when maximum Y is not in provided windows.
# TODO - [ ] Convenience function to quickly generate GRanges windows using syntax c("1:100-200", "X:500-800", "TP53", "1p31.1") with padding.

# Aesthetics
# TODO - [ ] SeqPlot-level aesthetic control of tracks and windows that defers to track-level changes to default.
# TODO - [ ] Separate axis lines, ticks, text, and titles into their own arguments (size, width, color, etc..)
# TODO - [ ] Aesthetic control for SeqElements (point/line size/color, etc)
# TODO - [ ] Y axis lines to windows
# TODO - [ ] Adding additional axis breaks (instead of just min and max)
# TODO - [ ] Adding panel grid breaks
# TODO - [ ] More margins to control spacing (between label and axis, between data and axis (expand))

# Bug Fixes
# TODO - [ ] When a track is plotted on top of a track containing SeqArches or SeqStrings, the y-scaling for the latter is altered.
