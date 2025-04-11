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
  upward <- orientation %in% c("*", "up", "+", "+/+", "+/-")
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

  upward <- orientation %in% c("*", "up", "+", "+/+", "+/-")
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
    windowBorderColor = "character"
  )
)

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
SeqTrack <- function(featureList, secondaryAxis = FALSE, xAxis = FALSE, yAxis = TRUE, yTitle = "", trackBackgroundColor = "transparent", trackBorderColor = "transparent", windowBackgroundColor = "transparent", windowBorderColor = "grey50") {
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
      windowBorderColor = windowBorderColor)
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
    for (feat in track@features) {
      if (inherits(feat, "SeqIdeogram")) {
        drawIdeogram(feat, layout_info, track_idx)
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
      if (is(sf, "SeqElement")) {
        yvals <- c(sf@y0, sf@y1)
        yvals[is.infinite(yvals)] <- NA
        if (class(sf) %in% c("SeqArch")) {
          yvals <- c(yvals, sf@height)
        }
        return(yvals)
      }
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
             cex = label_cex, srt = 90, adj = c(1, 0.5))
        text(x1, this_track_bottom - tick_length - 2, labels = paste0(round(scaled_end, 1), " ", unit_label),
             cex = label_cex, srt = 90, adj = c(1, 0.5))
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
        text(title_x, title_y, labels = track@yTitle, font = 2, adj = c(1, 0.5), cex = label_cex)
      }
    }
  }

  ## Process spanning features (SeqLink objects)
  for (track_idx in seq_len(totalTracks)) {
    track <- sp@tracks[[track_idx]]
    for (sf in track@features) {
      if (!inherits(sf, "SeqLink")) next
      # Use findOverlaps to map endpoints to windows
      ov1 <- findOverlaps(sf@gr1, global_windows)
      ov2 <- findOverlaps(sf@gr2, global_windows)
      valid_links <- intersect(queryHits(ov1), queryHits(ov2))
      if (length(valid_links) == 0) next

      for (i in valid_links) {
        w0 <- subjectHits(ov1)[queryHits(ov1) == i][1]
        w1 <- subjectHits(ov2)[queryHits(ov2) == i][1]
        if (is.na(w0) || is.na(w1)) next

        pt1 <- sf@gr1[i]
        pt2 <- sf@gr2[i]
        x0_abs <- globalTransformX(start(pt1), w0, layout_info)
        x1_abs <- globalTransformX(start(pt2), w1, layout_info)

        # Determine track indices for endpoints using t0/t1; default to current track if not set.
        track_index_0 <- if (!is.na(sf@t0[i]) && sf@t0[i] != 0) sf@t0[i] else track_idx
        track_index_1 <- if (!is.na(sf@t1[i]) && sf@t1[i] != 0) sf@t1[i] else track_idx

        # For arcs, if y0 and y1 are provided, use them; otherwise, default to the track's relative midpoint (0.5)
        y0_val <- sf@y0[i]
        y1_val <- sf@y1[i]
        if (!is.na(y0_val)) {
          y0_abs <- globalTransformY(y0_val, track_index_0, layout_info,
                                     track_y_min = rep(track_y_min, totalTracks),
                                     track_y_max = rep(track_y_max, totalTracks))
        } else {
          y0_abs <- (global_track_bottom <- layout_info$global_track_bottoms)[track_index_0] +
            0.5 * (layout_info$global_track_tops[track_index_0] - global_track_bottom[track_index_0])
        }
        if (!is.na(y1_val)) {
          y1_abs <- globalTransformY(y1_val, track_index_1, layout_info,
                                     track_y_min = rep(track_y_min, totalTracks),
                                     track_y_max = rep(track_y_max, totalTracks))
        } else {
          y1_abs <- (global_track_bottom <- layout_info$global_track_bottoms)[track_index_1] +
            0.5 * (layout_info$global_track_tops[track_index_1] - global_track_bottom[track_index_1])
        }

        # Retrieve the link's type and draw accordingly.
        if (inherits(sf, "SeqLine")) {
          segments(x0_abs, y0_abs, x1_abs, y1_abs, col = sf@color[i], lwd = 1.2)

        } else if (inherits(sf, "SeqArc")) {
          orientation <- if (!is.na(sf@orientation[i])) sf@orientation[i] else "*"
          #offset <- abs(x1_abs - x0_abs) * 0.25 * ifelse(orientation %in% c("*", "up", "+", "+/+", "+/-"), 1, -1)
          drawArc(x0_abs, y0_abs, x1_abs, y1_abs, orientation, col = sf@color[i], lwd = 1.2)

        } else if (is(sf, "SeqArch")) {
          arch_height <- sf@height[i]
          arch_orientation <- if (!is.na(sf@orientation[i])) sf@orientation[i] else "*"

          # Base of vertical lines (track-specific)
          track_index_0 <- if (!is.na(sf@t0[i]) && sf@t0[i] != 0) sf@t0[i] else track_idx
          track_index_1 <- if (!is.na(sf@t1[i]) && sf@t1[i] != 0) sf@t1[i] else track_idx

          y0_val <- sf@y0[i]
          y1_val <- sf@y1[i]

          y0_abs <- globalTransformY(y0_val, track_index_0, layout_info,
                                     track_y_min = rep(0, totalTracks),
                                     track_y_max = rep(1, totalTracks))

          y1_abs <- globalTransformY(y1_val, track_index_1, layout_info,
                                     track_y_min = rep(0, totalTracks),
                                     track_y_max = rep(1, totalTracks))

          # Top of vertical lines (original track) -
          arch_height <- sf@height[i]

          top0_abs <- globalTransformY(arch_height, track_idx, layout_info,
                                       track_y_min = rep(track_y_min, totalTracks),
                                       track_y_max = rep(track_y_max, totalTracks))

          top1_abs <- globalTransformY(arch_height, track_idx, layout_info,
                                       track_y_min = rep(track_y_min, totalTracks),
                                       track_y_max = rep(track_y_max, totalTracks))

          drawSeqArch(
            x0 = x0_abs,
            y0 = y0_abs,  # base of left stem (in t0)
            x1 = x1_abs,
            y1 = y1_abs,  # base of right stem (in t1)
            top0 = top0_abs,  # top of left stem (in current track)
            top1 = top1_abs,  # top of right stem (in current track)
            orientation = sf@orientation[i],
            curve = sf@curve[i],
            col = sf@color[i],
            lwd = 1.2
          )

        } else {
          warning("Unknown SeqLink subclass — skipping.")
        }

      }
    }
  }
})


# TODOS ----

# SeqFeatures
# TODO - [X] SeqFeature - Bar
# TODO - [ ] SeqFeature - Stacked Bar
# TODO - [ ] SeqFeature - Line (path through observations)
# TODO - [ ] SeqFeature - Area (basically line but a filled polygon)
# TODO - [ ] SeqFeature - Stacked Area
# TODO - [ ] SeqFeature - Segment
# TODO - [ ] SeqFeature - Tile (Heatmap with ranges on X and categorical on Y.)
# TODO - [ ] SeqFeature - Text

# SeqLinks (see genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1895-9 for nice examples)
# TODO - [X] SeqLink - Arc
# TODO - [X] SeqLink - Arch (optional height argument for y maximum, differs from arch in that it is bounded by vertical lines and added height functionality)
# TODO - [ ] SeqLink - ReCon (Special arch that separates types of rearrangements. Based on doi.org/10.1093/bioinformatics/btad719. Places inversions on one track and dup/del on another.
# TODO - [ ] SeqLink - String (optional height argument for y maximum)
# TODO - [ ] SeqLink - Barbell (intended to visualize paired reads)
# TODO - [ ] SeqLink - Band (for synteny-like connections)

# SeqAnnotations
# TODO - [X] SeqAnnotation - Ideogram
# TODO - [ ] SeqAnnotation - Gene/Transcript
# TODO - [ ] SeqAnnotation - Gene/Transcript
# TODO - [ ] SeqAnnotation - Box
# TODO - [ ] SeqAnnotation - Zoom
# TODO - [ ] SeqAnnotation - Text

# Window/Track Control
# TODO - [ ] Add functionality for customWindows at the Seq Track level. I.e., not all tracks have the same windows.
# TODO - [ ] Individual track and window gap control

# Aesthetics
# TODO - [ ] SeqPlot-level aesthetic control of tracks and windows that defers to track-level changes to default.
# TODO - [ ] Separate axis lines, ticks, text, and titles into their own arguments (size, width, color, etc..)
# TODO - [ ] Aesthetic control for SeqElements (point/line size/color, etc)
# TODO - [ ] Y axis lines to windows
# TODO - [ ] Adding additional axis breaks (instead of just min and max)
# TODO - [ ] Adding panel grid breaks
# TODO - [ ] More margins to control spacing (between label and axis, between data and axis (expand))
