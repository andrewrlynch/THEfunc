# SeqTile ----
#' SeqTile R6 Class
#'
#' @description
#' An R6 class for plotting categorical heatmap-style tiles.
#' Genomic ranges define the x-axis, and a categorical grouping column
#' defines the y-axis. Each tile is filled with the color stored in the
#' `color` metadata column of the GRanges.
#'
#' @export
SeqTile <- R6::R6Class("SeqTile",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic ranges.
                         gr = NULL,

                         #' @field groupCol Name of metadata column in `gr` used
                         #'   to define grouping (y-axis categories).
                         groupCol = NULL,

                         #' @field groups Vector of unique group labels, in display order.
                         groups = NULL,

                         #' @field y Integer vector mapping each row of `gr` to its
                         #'   group index (1..N). Used for discrete y positioning.
                         y = NULL,

                         #' @field coordOriginal Original GRanges object.
                         coordOriginal = NULL,

                         #' @field coordCanvas List of per-panel tile coordinates (x0, x1, y0, y1, color).
                         coordCanvas = NULL,

                         #' @field aesthetics Plotting aesthetics (border, lwd).
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default non-fill aesthetics for tiles.
                         defaultAesthetics = list(
                           border = NA,
                           lwd = 0.5
                         ),

                         #' @field style Visualization style: "full", "diagonal", "triangle", or "rectangle".
                         style = NULL,

                         #' @field yCoordType Y-axis coordinate type: "genomic" or "distance".
                         yCoordType = NULL,

                         #' @field maxDist Maximum genomic distance for clipping.
                         maxDist = NULL,

                         #' @description
                         #' Create a new `SeqTile` object.
                         #' Backwards compatibility: accepts `gr=` or `x=` for primary
                         #' genomic ranges. Optionally accepts `y=` as a `GRanges` to
                         #' enable genomic y-axis (2D mode). Otherwise use `aes(y=...)`
                         #' or `groupCol` to define discrete categories on the y-axis.
                         #' @param gr Deprecated primary `GRanges` (use `x=`).
                         #' @param x Primary `GRanges` for x-axis positions.
                         #' @param y Optional `GRanges` for genomic y-axis positions (must
                         #'   have length equal to `x` or `gr`).
                         #' @param groupCol Column name in `gr` metadata for grouping (y-axis).
                         #' @param aesthetics List of aesthetics (`border`, `lwd`).
                         #' @param aes A `SeqAes` object from `aes()`. Use `y` to map a
                         #'   column to discrete y-axis categories and `fill` to map a
                         #'   column to tile fill color.
                         #' @param scale A `SeqScale` color scale for the fill mapping.
                         #' @param style Visualization style: "full" (default), "diagonal",
                         #'   "triangle", or "rectangle". Only applies when `y=` GRanges
                         #'   is provided (2D genomic mode).
                         #' @param yCoordType Y-axis coordinate type: "genomic" (default) or
                         #'   "distance". Only used with rotated styles.
                         #' @param maxDist Maximum genomic distance for clipping. If NULL,
                         #'   defaults to max(width(y)) for "rectangle" style.
                         initialize = function(gr = NULL, x = NULL, y = NULL,
                                               groupCol = NULL, aesthetics = list(),
                                               aes = NULL, scale = NULL,
                                               style = "full", yCoordType = "genomic",
                                               maxDist = NULL) {
                           # resolve primary ranges (x or legacy gr)
                           if (!is.null(x)) {
                             stopifnot(inherits(x, "GRanges"))
                             gr_use <- x
                           } else if (!is.null(gr)) {
                             stopifnot(inherits(gr, "GRanges"))
                             gr_use <- gr
                           } else {
                             stop("Either `x=` (preferred) or `gr=` must be supplied as a GRanges for SeqTile.")
                           }

                           # --- Resolve y-axis grouping OR genomic y ---
                           if (!is.null(y)) {
                             # 2D genomic mode: y must be GRanges and match length
                             stopifnot(inherits(y, "GRanges"))
                             if (length(y) != length(gr_use))
                               stop("When providing `y=` GRanges it must be the same length as `x`/`gr`.")
                             self$gr_y <- y
                           } else if (!is.null(aes) && !is.null(aes[["y"]])) {
                             # discrete y mapping via aes
                             y_col <- aes[["y"]]
                             if (!y_col %in% names(mcols(gr_use)))
                               stop("y column '", y_col, "' not found in gr metadata.")
                             raw <- mcols(gr_use)[[y_col]]
                             self$groupCol <- y_col
                             self$groups   <- if (is.factor(raw)) levels(raw)
                                              else unique(as.character(raw))
                             self$y <- match(as.character(raw), self$groups)

                           } else if (!is.null(groupCol)) {
                             # Legacy path
                             if (!groupCol %in% names(mcols(gr_use)))
                               stop("Provided groupCol not found in metadata.")
                             self$groupCol <- groupCol
                             vals          <- as.character(mcols(gr_use)[[groupCol]])
                             self$groups   <- unique(vals)
                             self$y        <- match(vals, self$groups)
                           }

                           # --- Resolve fill ---
                           if (!is.null(aes) && !is.null(aes[["fill"]])) {
                             resolved <- .resolve_aes(
                               as.data.frame(S4Vectors::mcols(gr_use)), aes, scale,
                               length(gr_use), "grey80"
                             )
                             if (!is.null(resolved$fill)) {
                               mcols(gr_use)$color <- resolved$fill
                             }
                           } else if (!"color" %in% names(mcols(gr_use))) {
                             if (is.null(aes))
                               stop("`gr`/`x` must contain a 'color' metadata column or use aes(fill = ...).")
                           }

                           self$gr <- gr_use
                           self$coordOriginal <- gr_use
                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)

                           # Validate and store style parameters
                           style <- .validate_style_params(style, self$gr_y, maxDist, yCoordType)
                           self$style <- style
                           self$yCoordType <- yCoordType

                           # Auto-compute maxDist for rectangle style if needed
                           if (style == "rectangle" && is.null(maxDist)) {
                             if (!is.null(self$gr_y)) {
                               self$maxDist <- max(width(self$gr_y))
                             }
                           } else if (!is.null(maxDist)) {
                             self$maxDist <- maxDist
                           }
                         },

                         #' @description
                         #' Prepare tile coordinates for plotting in canvas space.
                         #' @param layout_track A list of panel metadata for the current track.
                         #' @param track_windows A `GRanges` of genomic windows for the track.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))

                           # For rectangle style, expand x-ranges before overlap detection
                           search_gr <- self$gr
                           if (self$style == "rectangle" && !is.null(self$gr_y) && !is.null(self$maxDist)) {
                             expanded_ranges <- ranges(self$gr)
                             expanded_ranges <- IRanges::IRanges(
                               start = pmax(start(expanded_ranges) - self$maxDist, 1),
                               end = end(expanded_ranges) + self$maxDist
                             )
                             search_gr <- GenomicRanges::GRanges(
                               seqnames = seqnames(self$gr),
                               ranges = expanded_ranges,
                               mcols = mcols(self$gr)
                             )
                           }

                           ov <- GenomicRanges::findOverlaps(search_gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- S4Vectors::queryHits(ov)
                           sh <- S4Vectors::subjectHits(ov)

                           # Apply diagonal filtering for 2D mode (before extracting data)
                           if (!is.null(self$gr_y) && self$style != "full") {
                             diag_mask <- .filter_diagonal_tiles(self$gr, self$gr_y, qh, sh, self$style)
                             qh <- qh[diag_mask]
                             sh <- sh[diag_mask]
                           }

                           # Extract data AFTER filtering (so indices match)
                           x0     <- start(self$gr)[qh]
                           x1     <- end(self$gr)[qh]
                           colors <- as.character(mcols(self$gr)$color)[qh]

                           # Store original coordinates before any transformations
                           x0_orig <- x0
                           x1_orig <- x1

                           # If genomic y provided, compute y0/y1 from gr_y
                           y_idx <- NULL
                           y0_raw <- NULL
                           y1_raw <- NULL
                           y0_orig <- NULL
                           y1_orig <- NULL
                           if (!is.null(self$gr_y)) {
                             y0_raw <- start(self$gr_y)[qh]
                             y1_raw <- end(self$gr_y)[qh]
                             y0_orig <- y0_raw
                             y1_orig <- y1_raw
                           } else if (!is.null(self$y)) {
                             y_idx <- self$y[qh]
                           }

                           # Pre-rotation clipping for rectangle style:
                           # Keep only tiles where x-range overlaps with ORIGINAL unexpanded windows
                           if (self$style == "rectangle" && !is.null(self$gr_y) && !is.null(self$maxDist)) {
                             # x0_orig, x1_orig are from unexpanded gr, not from expanded search_gr
                             # Check which tiles have x-range overlapping with track_windows
                             rect_mask <- rep(TRUE, length(qh))
                             for (i in seq_along(unique(sh))) {
                               w_idx <- unique(sh)[i]
                               window_range <- track_windows[w_idx]
                               w_mask <- sh == w_idx
                               # Keep tiles that overlap with actual window bounds
                               rect_mask[w_mask] <- end(self$gr)[qh[w_mask]] > start(window_range) &
                                                    start(self$gr)[qh[w_mask]] < end(window_range)
                             }
                             qh <- qh[rect_mask]
                             sh <- sh[rect_mask]
                             x0 <- x0[rect_mask]
                             x1 <- x1[rect_mask]
                             colors <- colors[rect_mask]
                             x0_orig <- x0_orig[rect_mask]
                             x1_orig <- x1_orig[rect_mask]
                             if (!is.null(y0_orig)) {
                               y0_orig <- y0_orig[rect_mask]
                               y1_orig <- y1_orig[rect_mask]
                             }
                             if (!is.null(y0_raw)) {
                               y0_raw <- y0_raw[rect_mask]
                               y1_raw <- y1_raw[rect_mask]
                             }
                             if (!is.null(y_idx)) {
                               y_idx <- y_idx[rect_mask]
                             }
                           }

                           for (w in unique(sh)) {
                             mask <- sh == w
                             if (sum(mask) == 0) next

                             panel <- layout_track[[w]]

                             # For rotated styles (triangle, rectangle), apply linear coordinate transformation
                             # in genomic space BEFORE canvas normalization: x_rot = (x+y)/2, y_rot = (y-x)/2
                             x0_work <- x0[mask]
                             x1_work <- x1[mask]

                             if (!is.null(y0_raw) && self$style %in% c("triangle", "rectangle")) {
                               # Apply linear coordinate transformation for 45° rotated styles
                               y0_work <- y0_raw[mask]
                               y1_work <- y1_raw[mask]

                               # Transform genomic (x, y) coordinates to rotated space
                               # x_rot = (x + y) / 2  (represents genomic position)
                               # y_rot = (y - x) / 2  (represents genomic distance)
                               x0_rot <- (x0_work + y0_work) / 2
                               x1_rot <- (x1_work + y1_work) / 2
                               y0_rot <- (y0_work - x0_work) / 2
                               y1_rot <- (y1_work - x1_work) / 2

                               x0_work <- x0_rot
                               x1_work <- x1_rot
                             } else {
                               y0_work <- NULL
                               y1_work <- NULL
                             }

                             # x: canvas mapping (works with original or transformed coordinates)
                             u0 <- (x0_work - panel$xscale[1]) / diff(panel$xscale)
                             u1 <- (x1_work - panel$xscale[1]) / diff(panel$xscale)
                             u0 <- pmax(pmin(u0, 1), 0)
                             u1 <- pmax(pmin(u1, 1), 0)

                             if (!is.null(y0_raw)) {
                               # Genomic y-axis
                               # For rotated styles, use already-transformed y0_work, y1_work
                               if (self$style %in% c("triangle", "rectangle")) {
                                 # For rotated styles, always use transformed distance coordinates
                                 # yCoordType is just metadata about axis representation (stored but not used for computation)
                                 y0_m <- y0_work
                                 y1_m <- y1_work
                               } else {
                                 # Full and diagonal styles: use original genomic coordinates
                                 y0_m <- y0_raw[mask]
                                 y1_m <- y1_raw[mask]
                               }

                               if (!is.null(panel$y_sub_panels)) {
                                 # Multiple y-windows: map tiles to appropriate sub-panel
                                 all_frames <- vector("list", length(panel$y_sub_panels))
                                 for (yw in seq_along(panel$y_sub_panels)) {
                                   sub <- panel$y_sub_panels[[yw]]
                                   yscale_yw <- sub$yscale
                                   y_in_sub <- y0_m >= yscale_yw[1] & y1_m <= yscale_yw[2]
                                   if (sum(y_in_sub) == 0) next
                                   v0 <- pmax(pmin((y0_m[y_in_sub] - yscale_yw[1]) / diff(yscale_yw), 1), 0)
                                   v1 <- pmax(pmin((y1_m[y_in_sub] - yscale_yw[1]) / diff(yscale_yw), 1), 0)
                                   all_frames[[yw]] <- data.frame(
                                     x0 = panel$inner$x0 + u0[y_in_sub] * (panel$inner$x1 - panel$inner$x0),
                                     x1 = panel$inner$x0 + u1[y_in_sub] * (panel$inner$x1 - panel$inner$x0),
                                     y0 = sub$y0 + v0 * (sub$y1 - sub$y0),
                                     y1 = sub$y0 + v1 * (sub$y1 - sub$y0),
                                     col = colors[mask][y_in_sub],
                                     x0_orig = x0_orig[mask][y_in_sub],
                                     x1_orig = x1_orig[mask][y_in_sub],
                                     y0_orig = y0_orig[mask][y_in_sub],
                                     y1_orig = y1_orig[mask][y_in_sub]
                                   )
                                 }
                                 self$coordCanvas[[w]] <- do.call(rbind, Filter(Negate(is.null), all_frames))
                               } else {
                                 # Single y-window
                                 v0 <- (y0_m - panel$yscale[1]) / diff(panel$yscale)
                                 v1 <- (y1_m - panel$yscale[1]) / diff(panel$yscale)
                                 v0 <- pmax(pmin(v0, 1), 0)
                                 v1 <- pmax(pmin(v1, 1), 0)

                                 # canvas coordinates
                                 x0_c <- panel$inner$x0 + u0 * (panel$inner$x1 - panel$inner$x0)
                                 x1_c <- panel$inner$x0 + u1 * (panel$inner$x1 - panel$inner$x0)
                                 y0_c <- panel$inner$y0 + v0 * (panel$inner$y1 - panel$inner$y0)
                                 y1_c <- panel$inner$y0 + v1 * (panel$inner$y1 - panel$inner$y0)

                                 self$coordCanvas[[w]] <- data.frame(
                                   x0 = x0_c, x1 = x1_c,
                                   y0 = y0_c, y1 = y1_c,
                                   col = colors[mask],
                                   x0_orig = x0_orig[mask],
                                   x1_orig = x1_orig[mask],
                                   y0_orig = y0_orig[mask],
                                   y1_orig = y1_orig[mask]
                                 )
                               }
                             } else {
                               # y: integer positions scaled through yscale
                               # For discrete yscale = c(0.5, N+0.5),
                               # each tile spans (idx - 0.5) to (idx + 0.5)
                               y0_data <- y_idx[mask] - 0.5
                               y1_data <- y_idx[mask] + 0.5
                               v0 <- (y0_data - panel$yscale[1]) / diff(panel$yscale)
                               v1 <- (y1_data - panel$yscale[1]) / diff(panel$yscale)
                               v0 <- pmax(pmin(v0, 1), 0)
                               v1 <- pmax(pmin(v1, 1), 0)

                               # canvas coordinates
                               x0_c <- panel$inner$x0 + u0 * (panel$inner$x1 - panel$inner$x0)
                               x1_c <- panel$inner$x0 + u1 * (panel$inner$x1 - panel$inner$x0)
                               y0_c <- panel$inner$y0 + v0 * (panel$inner$y1 - panel$inner$y0)
                               y1_c <- panel$inner$y0 + v1 * (panel$inner$y1 - panel$inner$y0)

                               self$coordCanvas[[w]] <- data.frame(
                                 x0 = x0_c, x1 = x1_c,
                                 y0 = y0_c, y1 = y1_c,
                                 col = colors[mask],
                                 x0_orig = x0_orig[mask],
                                 x1_orig = x1_orig[mask],
                                 y0_orig = rep(NA_real_, sum(mask)),
                                 y1_orig = rep(NA_real_, sum(mask))
                               )
                             }
                           }

                           invisible()
                         },

                         #' @description
                         #' Draw the tiles using grid graphics.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()
                           for (w in seq_along(self$coordCanvas)) {
                             coords <- self$coordCanvas[[w]]
                             if (is.null(coords) || nrow(coords) == 0) next

                             if (self$style %in% c("full", "diagonal")) {
                               # Use grid.rect() for full and diagonal styles
                               grid.rect(
                                 x = unit((coords$x0 + coords$x1) / 2, "npc"),
                                 y = unit((coords$y0 + coords$y1) / 2, "npc"),
                                 width  = unit(coords$x1 - coords$x0, "npc"),
                                 height = unit(coords$y1 - coords$y0, "npc"),
                                 gp = gpar(
                                   fill = coords$col,
                                   col = self$aesthetics$border,
                                   lwd = self$aesthetics$lwd
                                 ),
                                 just = c("center", "center")
                               )
                             } else if (self$style %in% c("triangle", "rectangle")) {
                               # For rotated styles (triangle, rectangle), coordinates in coordCanvas
                               # are already in rotated space (after linear transformation in prep()).
                               # We render them as diamonds using grid.polygon() with 4 corners.
                               # The corners form a diamond shape in the rotated coordinate space.
                               for (i in seq_len(nrow(coords))) {
                                 # Compute 4 corners of the diamond in canvas space
                                 # The rectangle (x0, y0) to (x1, y1) in rotated space
                                 # becomes a diamond with 4 corners
                                 x_center <- (coords$x0[i] + coords$x1[i]) / 2
                                 y_center <- (coords$y0[i] + coords$y1[i]) / 2
                                 x_half <- (coords$x1[i] - coords$x0[i]) / 2
                                 y_half <- (coords$y1[i] - coords$y0[i]) / 2

                                 # Diamond corners (left, bottom, right, top)
                                 diamond_x <- c(x_center - x_half, x_center, x_center + x_half, x_center)
                                 diamond_y <- c(y_center, y_center - y_half, y_center, y_center + y_half)

                                 grid.polygon(
                                   x = unit(diamond_x, "npc"),
                                   y = unit(diamond_y, "npc"),
                                   gp = gpar(
                                     fill = coords$col[i],
                                     col = self$aesthetics$border,
                                     lwd = self$aesthetics$lwd
                                   )
                                 )
                               }
                             }
                           }
                         },

                         .infer_scale_y = function() {
                           if (!is.null(self$groups))
                             seq_scale_discrete(levels = self$groups)
                           else
                             NULL
                         }
                       )
)

# Helper functions for SeqTile style support

#' @keywords internal
.validate_style_params <- function(style, gr_y, maxDist, yCoordType) {
  valid_styles <- c("full", "diagonal", "triangle", "rectangle")
  if (!style %in% valid_styles) {
    stop("style must be one of: ", paste(valid_styles, collapse = ", "))
  }

  if (style != "full" && is.null(gr_y)) {
    warning("style='", style, "' applies only with 2D genomic y-axis (y=GRanges). Reverting to 'full'.")
    return("full")
  }

  if (!is.null(maxDist)) {
    if (!is.numeric(maxDist) || maxDist <= 0) {
      stop("maxDist must be a positive number")
    }
  }

  if (!yCoordType %in% c("genomic", "distance")) {
    stop("yCoordType must be 'genomic' or 'distance'")
  }

  style
}

#' @keywords internal
.filter_diagonal_tiles <- function(gr_x, gr_y, qh, sh, style) {
  if (style == "full") {
    return(rep(TRUE, length(qh)))
  }

  if (style == "diagonal") {
    # Keep upper diagonal: y >= x (in genomic coordinates)
    y_start <- start(gr_y)[qh]
    x_start <- start(gr_x)[qh]
    return(y_start >= x_start)
  }

  if (style %in% c("triangle", "rectangle")) {
    # Keep upper-right triangle region
    y_end <- end(gr_y)[qh]
    x_start <- start(gr_x)[qh]
    y_start <- start(gr_y)[qh]
    x_end <- end(gr_x)[qh]
    return(y_end >= x_start & y_start <= x_end)
  }

  rep(TRUE, length(qh))
}

#' @keywords internal
.transform_to_rotated_coords <- function(x0, x1, y0, y1) {
  # Linear coordinate transformation for 45° rotated Hi-C style heatmaps
  # Based on the approach from hic_rotated_heatmap.R
  #
  # Transformation:
  #   x_rot = (x + y) / 2  (genomic position - bin average)
  #   y_rot = (y - x) / 2  (genomic distance - half bin difference)
  #
  # For a rectangular tile with corners (x0,y0), (x1,y0), (x1,y1), (x0,y1),
  # the transformed corners form a diamond shape in rotated space.
  # Corner transformation:
  #   (x0, y0) → ((x0+y0)/2, (y0-x0)/2)  - left
  #   (x1, y0) → ((x1+y0)/2, (y0-x1)/2)  - bottom
  #   (x1, y1) → ((x1+y1)/2, (y1-x1)/2)  - right
  #   (x0, y1) → ((x0+y1)/2, (y1-x0)/2)  - top

  list(
    x = c(
      (x0 + y0) / 2,  # left
      (x1 + y0) / 2,  # bottom
      (x1 + y1) / 2,  # right
      (x0 + y1) / 2   # top
    ),
    y = c(
      (y0 - x0) / 2,  # left
      (y0 - x1) / 2,  # bottom
      (y1 - x1) / 2,  # right
      (y1 - x0) / 2   # top
    )
  )
}