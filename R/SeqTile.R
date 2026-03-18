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

                         #' @field panelBounds List of per-panel inner-panel NPC bounds used for
                         #'   viewport clipping in rotated styles (triangle, rectangle).
                         panelBounds = NULL,

                         #' @field aesthetics Plotting aesthetics (border, lwd).
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default non-fill aesthetics for tiles.
                         defaultAesthetics = list(
                           border = NA,
                           lwd = 0.5
                         ),

                         #' @field style Visualization style: "full", "diagonal", "triangle", or "rectangle".
                         style = NULL,

                         #' @field yCoordType Y-axis coordinate type: auto-set to "distance" for triangle/rectangle styles, "genomic" for full/diagonal.
                         yCoordType = NULL,

                         #' @field maxDist Maximum genomic distance for clipping.
                         maxDist = NULL,

                         #' @field yDistMax Maximum distance (bp) for the y-axis in triangle/rectangle styles.
                         #'   If NULL, falls back to maxDist. Only used with rotated styles.
                         yDistMax = NULL,

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
                         #' @param maxDist Maximum genomic distance for clipping. If NULL,
                         #'   defaults to max(width(y)) for "rectangle" style.
                         #' @param yDistMax Maximum distance (bp) shown on the y-axis for
                         #'   triangle/rectangle styles. If NULL, falls back to `maxDist`.
                         initialize = function(gr = NULL, x = NULL, y = NULL,
                                               groupCol = NULL, aesthetics = list(),
                                               aes = NULL, scale = NULL,
                                               style = "full",
                                               maxDist = NULL, yDistMax = NULL) {
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
                           style <- .validate_style_params(style, self$gr_y, maxDist)
                           self$style <- style
                           self$yCoordType <- if (style %in% c("triangle", "rectangle")) "distance" else "genomic"

                           # Auto-compute maxDist for rectangle style if needed
                           if (style == "rectangle" && is.null(maxDist)) {
                             if (!is.null(self$gr_y)) {
                               self$maxDist <- max(width(self$gr_y))
                             }
                           } else if (!is.null(maxDist)) {
                             self$maxDist <- maxDist
                           }
                           self$yDistMax <- yDistMax
                         },

                         #' @description
                         #' Prepare tile coordinates for plotting in canvas space.
                         #' @param layout_track A list of panel metadata for the current track.
                         #' @param track_windows A `GRanges` of genomic windows for the track.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))
                           self$panelBounds <- vector("list", length(track_windows))

                           # For rectangle style, expand x-ranges before overlap detection.
                           # Tiles at the upper-left corner have x_orig ≈ window_min - yDistMax/2,
                           # which is OUTSIDE the x-window. Expanding the search range by the
                           # max visible distance ensures these tiles are loaded so polygon
                           # clipping in draw() can fill the full upper-left corner.
                           # Use maxDist if set, else yDistMax, else the largest window width.
                           search_gr <- self$gr
                           if (self$style == "rectangle" && !is.null(self$gr_y)) {
                             expand_by <- self$maxDist %||% self$yDistMax %||%
                                          max(GenomicRanges::width(track_windows))
                             expanded_ranges <- IRanges::IRanges(
                               start = pmax(start(ranges(self$gr)) - expand_by, 1),
                               end   = end(ranges(self$gr)) + expand_by
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

                           # Rectangle: also search by y-axis overlap to capture upper-left corner tiles.
                           # Those tiles have x_orig < window_min (absent from x-based findOverlaps when
                           # Hi-C data was loaded only for the visible x-window) but y_orig inside the
                           # window. The upper-right corner is already covered because those tiles have
                           # x_orig inside the window. The straddle filter will discard any tiles that
                           # fall entirely outside the x_rot range, so no extra data leaks in.
                           if (self$style == "rectangle" && !is.null(self$gr_y)) {
                             ov_y <- GenomicRanges::findOverlaps(self$gr_y, track_windows)
                             if (length(ov_y) > 0) {
                               qh_all <- c(qh, S4Vectors::queryHits(ov_y))
                               sh_all <- c(sh, S4Vectors::subjectHits(ov_y))
                               keep   <- !duplicated(paste(qh_all, sh_all, sep = "_"))
                               qh     <- qh_all[keep]
                               sh     <- sh_all[keep]
                             }
                           }

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

                           # NOTE: For "rectangle" style, do NOT apply any pre-rotation clipping here.
                           # The window expansion above already fetched flanking contacts (x outside
                           # the display window but whose rotated position x_rot = (x+y)/2 lands
                           # inside the window). Clipping them back to the original x-range would
                           # remove exactly the contacts that fill the top corners of the rectangle.
                           # The u0/u1 clamping in the coordinate normalisation below handles the
                           # visual boundary cleanly.

                           for (w in unique(sh)) {
                             mask <- sh == w
                             if (sum(mask) == 0) next

                             panel <- layout_track[[w]]

                            # Filter out tiles outside the diagonal edges (for triangle/rectangle styles).
                            # CRITICAL: Use UNEXPANDED bounds (panel$data_x) for clipping decisions,
                            # not expanded bounds (panel$xscale). Expansion is for viewport display only.
                            # NOTE: left_ok/right_ok have length sum(mask), not length(mask),
                            # so must use mask[mask] <- ... to avoid silent R vector recycling.
                            if (self$style %in% c("triangle", "rectangle") && !is.null(y0_orig)) {
                              data_x <- panel$data_x %||% panel$xscale  # fallback if data_x not set

                              if (self$style == "triangle") {
                                # Triangle: tile omission creates clean diagonal edges.
                                # Simple y-coord check works because Hi-C is symmetric (same chr, same scale).
                                # Removes ANY tile that has y < window_min or y > window_max.
                                left_ok  <- !(y0_orig[mask] < data_x[1])
                                right_ok <- !(y1_orig[mask] > data_x[2])
                              } else {
                                # Rectangle: straddling tiles at the diagonal boundaries must be KEPT
                                # so polygon clipping in draw() can produce clean diagonal edges.
                                # Only remove tiles with NO x_rot overlap with [window_min, window_max]:
                                #   x_rot range of a tile = [(x0+y0)/2, (x1+y1)/2]
                                #   Entirely left  → (x1+y1)/2 < window_min → x1+y1 < 2*window_min
                                #   Entirely right → (x0+y0)/2 > window_max → x0+y0 > 2*window_max
                                left_ok  <- !((x1_orig[mask] + y1_orig[mask]) < 2 * data_x[1])
                                right_ok <- !((x0_orig[mask] + y0_orig[mask]) > 2 * data_x[2])
                              }

                              mask[mask] <- left_ok & right_ok
                              if (sum(mask) == 0) next
                            }

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
                               #
                               # The four corners of a rectangular tile (x0,y0)-(x1,y1) become
                               # a diamond with corners:
                               #   left:   ((x0+y0)/2, (y0-x0)/2)
                               #   bottom: ((x1+y0)/2, (y0-x1)/2)
                               #   right:  ((x1+y1)/2, (y1-x1)/2)
                               #   top:    ((x0+y1)/2, (y1-x0)/2)
                               # The bounding-box x range  = [left.x,  right.x]  = [(x0+y0)/2, (x1+y1)/2]
                               # The bounding-box y range  = [bottom.y, top.y]   = [(y0-x1)/2, (y1-x0)/2]
                               # (uses CROSS terms for y so that square tiles have non-zero height)
                               x0_rot <- (x0_work + y0_work) / 2
                               x1_rot <- (x1_work + y1_work) / 2
                               y0_rot <- (y0_work - x1_work) / 2  # bottom corner y
                               y1_rot <- (y1_work - x0_work) / 2  # top corner y

                               x0_work <- x0_rot
                               x1_work <- x1_rot
                               y0_work <- y0_rot
                               y1_work <- y1_rot
                             } else {
                               y0_work <- NULL
                               y1_work <- NULL
                             }

                             # x: canvas mapping (works with original or transformed coordinates)
                             u0 <- (x0_work - panel$xscale[1]) / diff(panel$xscale)
                             u1 <- (x1_work - panel$xscale[1]) / diff(panel$xscale)
                             u0_raw <- u0
                             u1_raw <- u1
                             # Rotated styles need unclamped u for clip detection; clamp others now.
                             # Clamp to DATA boundaries (not expanded xscale) so tiles don't bleed
                             # into the xExpansion gap. Fall back to [0,1] when data_x is absent.
                             if (!self$style %in% c("triangle", "rectangle")) {
                               if (!is.null(panel$data_x) && diff(panel$xscale) > 0) {
                                 u_lo <- (panel$data_x[1] - panel$xscale[1]) / diff(panel$xscale)
                                 u_hi <- (panel$data_x[2] - panel$xscale[1]) / diff(panel$xscale)
                               } else {
                                 u_lo <- 0; u_hi <- 1
                               }
                               u0 <- pmax(pmin(u0, u_hi), u_lo)
                               u1 <- pmax(pmin(u1, u_hi), u_lo)
                             }

                             if (!is.null(y0_raw)) {
                               # Genomic y-axis
                               # For rotated styles, use already-transformed y0_work, y1_work
                               if (self$style %in% c("triangle", "rectangle")) {
                                 # For rotated styles, use transformed distance coordinates (y_rot)
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
                                 # Single y-window — determine effective yscale and v0/v1.

                                 if (self$style %in% c("triangle", "rectangle")) {
                                   # Distance mode: use FULL distance (2 × half-distance y_rot).
                                   # Aligns with the [0, yDistMax] scale from .infer_scale_y().
                                   y0_m <- y0_m * 2   # = y0_raw - x1_raw (full distance, bp)
                                   y1_m <- y1_m * 2   # = y1_raw - x0_raw (full distance, bp)
                                   dist_max <- self$yDistMax %||% self$maxDist %||%
                                               diff(panel$data_x %||% panel$xscale)

                                   # Apply yExpansion proportionally to the distance scale.
                                   # Derive expansion fractions from the panel's expanded vs unexpanded
                                   # yscale bounds, then apply proportionally to dist_max.
                                   # This correctly handles BOTH cases:
                                   #   - yDistMax set:     panel$data_y = c(0, dist_max), fracs are correct
                                   #   - yDistMax not set: panel$data_y = c(0, window_width) (set by
                                   #     SeqPlot's new fallback), fracs are also correct
                                   if (!is.null(panel$yscale) && !is.null(panel$data_y)) {
                                     span_panel <- diff(panel$data_y)
                                     if (span_panel > 0) {
                                       frac_lower <- max(0, panel$data_y[1] - panel$yscale[1]) / span_panel
                                       frac_upper <- max(0, panel$yscale[2] - panel$data_y[2]) / span_panel
                                       lower_pad <- frac_lower * dist_max
                                       upper_pad <- frac_upper * dist_max
                                     } else {
                                       lower_pad <- 0
                                       upper_pad <- 0
                                     }
                                     yscale_eff <- c(0 - lower_pad, dist_max + upper_pad)
                                   } else {
                                     # No expansion info available, use plain distance scale
                                     yscale_eff <- c(0, dist_max)
                                   }
                                 } else {
                                   # Genomic mode for full/diagonal styles.
                                   # Use panel$yscale when a meaningful y scale was provided
                                   # (y_windows or infer_scale_y()). panel$yscale already
                                   # includes yExpansion from .expand_limits() in SeqPlot.
                                   # Tiles with y outside yscale_eff are clamped to
                                   # v0=v1=boundary (zero height → invisible).
                                   # Fall back to data range only when panel$data_y is the
                                   # default c(0,1) sentinel, meaning no explicit y scale
                                   # was set.
                                   yscale_eff <- if (!is.null(panel$data_y) &&
                                                     !(panel$data_y[1] == 0 &&
                                                       panel$data_y[2] == 1)) {
                                     panel$yscale   # y_windows / infer_scale_y → includes yExpansion
                                   } else {
                                     range(c(y0_m, y1_m))   # no explicit scale → auto-range
                                   }
                                 }

                                 v0     <- (y0_m - yscale_eff[1]) / diff(yscale_eff)
                                 v1     <- (y1_m - yscale_eff[1]) / diff(yscale_eff)
                                 v0_raw <- v0
                                 v1_raw <- v1
                                 # Clamp for non-rotated styles only.
                                 # Clamp to DATA boundaries (not expanded yscale) so tiles don't
                                 # bleed into the yExpansion gap.  Tiles entirely outside the data
                                 # range collapse to zero height and become invisible.  Tiles
                                 # straddling the boundary are clipped to start/end at the boundary.
                                 if (!self$style %in% c("triangle", "rectangle")) {
                                   if (!is.null(panel$data_y) &&
                                       !(panel$data_y[1] == 0 && panel$data_y[2] == 1) &&
                                       diff(yscale_eff) > 0) {
                                     v_lo <- (panel$data_y[1] - yscale_eff[1]) / diff(yscale_eff)
                                     v_hi <- (panel$data_y[2] - yscale_eff[1]) / diff(yscale_eff)
                                   } else {
                                     v_lo <- 0; v_hi <- 1
                                   }
                                   v0 <- pmax(pmin(v0, v_hi), v_lo)
                                   v1 <- pmax(pmin(v1, v_hi), v_lo)
                                 }

                                 # Canvas coordinates — keep unclamped so the bounding box
                                 # reflects the true diamond shape for correct clip geometry.
                                 x0_c <- panel$inner$x0 + u0 * (panel$inner$x1 - panel$inner$x0)
                                 x1_c <- panel$inner$x0 + u1 * (panel$inner$x1 - panel$inner$x0)
                                 y0_c <- panel$inner$y0 + v0 * (panel$inner$y1 - panel$inner$y0)
                                 y1_c <- panel$inner$y0 + v1 * (panel$inner$y1 - panel$inner$y0)

                                 if (self$style %in% c("triangle", "rectangle")) {
                                   # Clip-type bit flags: 1=left, 2=right, 4=bottom, 8=top
                                   clip <- bitwOr(
                                     bitwOr(as.integer(u0_raw < 0) * 1L,
                                            as.integer(u1_raw > 1) * 2L),
                                     bitwOr(as.integer(v0_raw < 0) * 4L,
                                            as.integer(v1_raw > 1) * 8L)
                                   )
                                   # Discard tiles fully outside the visible area
                                   visible <- u1_raw > 0 & u0_raw < 1 &
                                              v1_raw > 0 & v0_raw < 1

                                   # Store inner-panel NPC bounds for viewport clipping in draw().
                                   ix0 <- panel$inner$x0; ix1 <- panel$inner$x1
                                   iy0 <- panel$inner$y0; iy1 <- panel$inner$y1
                                   self$panelBounds[[w]] <- list(
                                     x0 = ix0, x1 = ix1,
                                     y0 = iy0, y1 = iy1,
                                     xscale = panel$xscale,
                                     yscale = yscale_eff,
                                     # Unexpanded data range in distance coords [0, dist_max].
                                     # Used in draw() to clip at data boundaries (not physical
                                     # panel edges) so yExpansion creates proper gaps.
                                     data_yscale = if (self$style %in% c("triangle", "rectangle"))
                                       c(0, dist_max) else NULL,
                                     # Unexpanded x window range [window_min, window_max].
                                     # Used in rectangle draw() to compute diagonal clip NPC positions.
                                     data_xscale = panel$data_x
                                   )

                                   self$coordCanvas[[w]] <- data.frame(
                                     x0  = x0_c[visible], x1 = x1_c[visible],
                                     y0  = y0_c[visible], y1 = y1_c[visible],
                                     col = colors[mask][visible],
                                     x0_orig = x0_orig[mask][visible],
                                     x1_orig = x1_orig[mask][visible],
                                     y0_orig = y0_orig[mask][visible],
                                     y1_orig = y1_orig[mask][visible]
                                   )
                                 } else {
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
                               x0  <- coords$x0
                               x1  <- coords$x1
                               y0  <- coords$y0
                               y1  <- coords$y1
                               col <- coords$col
                               brd <- self$aesthetics$border
                               lwd <- self$aesthetics$lwd

                               pb  <- if (!is.null(self$panelBounds)) self$panelBounds[[w]] else NULL
                               has_bounds <- !is.null(pb) && (pb$x1 > pb$x0) && (pb$y1 > pb$y0)

                               xc <- (x0 + x1) / 2  # unclamped diamond centres
                               yc <- (y0 + y1) / 2

                               if (has_bounds && self$style == "rectangle") {
                                 # Rectangle: clip each diamond against all 4 data boundaries.
                                 # Clip at DATA boundaries (not physical panel edges) so:
                                 #   - yExpansion creates proper gaps at top/bottom
                                 #   - xExpansion doesn't allow tiles to bleed into x-margins
                                 #   - Diagonal edges (left/right in rotated space) come from
                                 #     polygon clipping, not tile omission (unlike triangle).
                                 span_npc_x <- pb$x1 - pb$x0
                                 span_npc_y <- pb$y1 - pb$y0
                                 x_span <- diff(pb$xscale)
                                 y_span <- diff(pb$yscale)

                                 # X clip: window_min and window_max → the diagonal edges
                                 if (!is.null(pb$data_xscale) && x_span > 0) {
                                   x_left  <- pb$x0 + (pb$data_xscale[1] - pb$xscale[1]) / x_span * span_npc_x
                                   x_right <- pb$x0 + (pb$data_xscale[2] - pb$xscale[1]) / x_span * span_npc_x
                                 } else {
                                   x_left <- pb$x0; x_right <- pb$x1
                                 }

                                 # Y clip: distance=0 (bottom) and dist_max (top) → horizontal edges
                                 if (!is.null(pb$data_yscale) && y_span > 0) {
                                   y_bottom <- pb$y0 + (pb$data_yscale[1] - pb$yscale[1]) / y_span * span_npc_y
                                   y_top    <- pb$y0 + (pb$data_yscale[2] - pb$yscale[1]) / y_span * span_npc_y
                                 } else {
                                   y_bottom <- pb$y0; y_top <- pb$y1
                                 }

                                 for (i in seq_along(x0)) {
                                   diamond_x <- c(x0[i], xc[i], x1[i], xc[i], x0[i])  # closed polygon
                                   diamond_y <- c(yc[i], y0[i], yc[i], y1[i], yc[i])

                                   clipped <- self$.clip_polygon_rect(
                                     diamond_x, diamond_y,
                                     xmin = x_left, xmax = x_right,
                                     ymin = y_bottom, ymax = y_top
                                   )

                                   if (length(clipped$x) > 0) {
                                     grid.polygon(
                                       x  = unit(clipped$x, "npc"),
                                       y  = unit(clipped$y, "npc"),
                                       gp = gpar(fill = col[i], col = brd, lwd = lwd)
                                     )
                                   }
                                 }
                               } else if (has_bounds && self$style == "triangle") {
                                 # Triangle: tile omission (straddling filter in prep()) handles diagonal
                                 # left/right edges. Horizontal polygon clipping handles top/bottom edges.
                                 #
                                 # Clip at DATA boundaries (distance=0, distance=dist_max), NOT at
                                 # physical panel edges (pb$y0, pb$y1). This ensures yExpansion creates
                                 # proper gaps: distance=0 sits above pb$y0, dist_max sits below pb$y1.
                                 span_npc <- pb$y1 - pb$y0
                                 dist_span <- diff(pb$yscale)
                                 if (!is.null(pb$data_yscale) && dist_span > 0) {
                                   y_bottom <- pb$y0 + (pb$data_yscale[1] - pb$yscale[1]) / dist_span * span_npc
                                   y_top    <- pb$y0 + (pb$data_yscale[2] - pb$yscale[1]) / dist_span * span_npc
                                 } else {
                                   y_bottom <- pb$y0  # fallback (no expansion info)
                                   y_top    <- pb$y1
                                 }
                                 for (i in seq_along(x0)) {
                                   diamond_x <- c(x0[i], xc[i], x1[i], xc[i])
                                   diamond_y <- c(yc[i], y0[i], yc[i], y1[i])
                                   # Clip at top (y = y_top = dist_max in NPC)
                                   clipped <- self$.clip_polygon_horizontal_edge(diamond_x, diamond_y, y_top, keep_above = FALSE)
                                   # Clip at bottom (y = y_bottom = distance 0 in NPC)
                                   if (length(clipped$x) > 0)
                                     clipped <- self$.clip_polygon_horizontal_edge(clipped$x, clipped$y, y_bottom, keep_above = TRUE)
                                   if (length(clipped$x) > 0) {
                                     grid.polygon(
                                       x  = unit(clipped$x, "npc"),
                                       y  = unit(clipped$y, "npc"),
                                       gp = gpar(fill = col[i], col = brd, lwd = lwd)
                                     )
                                   }
                                 }
                               } else {
                                 # No clipping needed: draw all diamonds as-is
                                 grid.polygon(
                                   x          = unit(c(rbind(x0, xc, x1, xc)), "npc"),
                                   y          = unit(c(rbind(yc, y0, yc, y1)), "npc"),
                                   id.lengths = rep(4L, length(x0)),
                                   gp         = gpar(fill = col, col = brd, lwd = lwd)
                                 )
                               }
                             }
                           }
                         },

                         # Clipping helpers for polygon-based viewport clipping

                        #' @description
                        #' Clip a convex polygon against a rectangular clipping region using
                        #' Sutherland-Hodgman algorithm. Returns a list of (x, y) vertex vectors.
                        #' Used for "rectangle" style to create straight-edged clipping.
                        #'
                        #' @param x,y  Numeric vectors of polygon vertices (must be closed)
                        #' @param xmin,xmax,ymin,ymax  Rectangular clip region bounds
                        #'
                        #' @return List(x, y) with clipped polygon vertices, possibly fewer
                        .clip_polygon_rect = function(x, y, xmin, xmax, ymin, ymax) {
                          # Sutherland-Hodgman: clip polygon against 4 rectangular edges in sequence
                          clip_left   = function(x, y) self$.clip_polygon_edge(x, y, "left",   xmin)
                          clip_right  = function(x, y) self$.clip_polygon_edge(x, y, "right",  xmax)
                          clip_bottom = function(x, y) self$.clip_polygon_edge(x, y, "bottom", ymin)
                          clip_top    = function(x, y) self$.clip_polygon_edge(x, y, "top",    ymax)

                          result <- list(x = x, y = y)
                          for (clip_fn in list(clip_left, clip_right, clip_bottom, clip_top)) {
                            result <- clip_fn(result$x, result$y)
                            if (length(result$x) == 0) break
                          }
                          result
                        },

                        #' @description
                        #' Clip polygon against a single axis-aligned edge.
                        .clip_polygon_edge = function(x, y, edge, val) {
                          if (length(x) == 0) return(list(x = numeric(0), y = numeric(0)))

                          if (x[1] != x[length(x)] || y[1] != y[length(y)]) {
                            x <- c(x, x[1])
                            y <- c(y, y[1])
                          }

                          n <- length(x) - 1
                          out_x <- numeric(0)
                          out_y <- numeric(0)

                          if (edge == "left") {
                            inside <- function(xx) xx >= val
                            intersect <- function(x1, y1, x2, y2) {
                              if (x1 == x2) return(list(x = x1, y = y1))
                              t <- (val - x1) / (x2 - x1)
                              list(x = val, y = y1 + t * (y2 - y1))
                            }
                          } else if (edge == "right") {
                            inside <- function(xx) xx <= val
                            intersect <- function(x1, y1, x2, y2) {
                              if (x1 == x2) return(list(x = x1, y = y1))
                              t <- (val - x1) / (x2 - x1)
                              list(x = val, y = y1 + t * (y2 - y1))
                            }
                          } else if (edge == "bottom") {
                            inside <- function(yy) yy >= val
                            intersect <- function(x1, y1, x2, y2) {
                              if (y1 == y2) return(list(x = x1, y = y1))
                              t <- (val - y1) / (y2 - y1)
                              list(x = x1 + t * (x2 - x1), y = val)
                            }
                          } else if (edge == "top") {
                            inside <- function(yy) yy <= val
                            intersect <- function(x1, y1, x2, y2) {
                              if (y1 == y2) return(list(x = x1, y = y1))
                              t <- (val - y1) / (y2 - y1)
                              list(x = x1 + t * (x2 - x1), y = val)
                            }
                          }

                          for (i in 1:n) {
                            x1 <- x[i];  y1 <- y[i]
                            x2 <- x[i+1]; y2 <- y[i+1]
                            inside1 <- if (edge %in% c("left", "right")) inside(x1) else inside(y1)
                            inside2 <- if (edge %in% c("left", "right")) inside(x2) else inside(y2)

                            if (inside2) {
                              if (!inside1) {
                                inter <- intersect(x1, y1, x2, y2)
                                out_x <- c(out_x, inter$x)
                                out_y <- c(out_y, inter$y)
                              }
                              out_x <- c(out_x, x2)
                              out_y <- c(out_y, y2)
                            } else if (inside1) {
                              inter <- intersect(x1, y1, x2, y2)
                              out_x <- c(out_x, inter$x)
                              out_y <- c(out_y, inter$y)
                            }
                          }

                          list(x = out_x, y = out_y)
                        },

                        #' @description
                        #' Clip polygon against a horizontal line y = y_val.
                        #' keep_above=TRUE keeps region y >= y_val, FALSE keeps y <= y_val.
                        .clip_polygon_horizontal_edge = function(x, y, y_val, keep_above = FALSE) {
                          if (length(x) == 0) return(list(x = numeric(0), y = numeric(0)))

                          if (x[1] != x[length(x)] || y[1] != y[length(y)]) {
                            x <- c(x, x[1])
                            y <- c(y, y[1])
                          }

                          n <- length(x) - 1
                          out_x <- numeric(0)
                          out_y <- numeric(0)

                          inside <- function(yi) {
                            if (keep_above) yi >= y_val else yi <= y_val
                          }

                          for (i in 1:n) {
                            x1 <- x[i];   y1 <- y[i]
                            x2 <- x[i+1]; y2 <- y[i+1]
                            i1 <- inside(y1)
                            i2 <- inside(y2)

                            if (i2) {
                              if (!i1) {
                                # Crossing from outside to inside: interpolate
                                if (abs(y2 - y1) > 1e-12) {
                                  t <- (y_val - y1) / (y2 - y1)
                                  out_x <- c(out_x, x1 + t * (x2 - x1))
                                } else {
                                  out_x <- c(out_x, x1)
                                }
                                out_y <- c(out_y, y_val)
                              }
                              out_x <- c(out_x, x2)
                              out_y <- c(out_y, y2)
                            } else if (i1) {
                              # Crossing from inside to outside: interpolate
                              if (abs(y2 - y1) > 1e-12) {
                                t <- (y_val - y1) / (y2 - y1)
                                out_x <- c(out_x, x1 + t * (x2 - x1))
                              } else {
                                out_x <- c(out_x, x1)
                              }
                              out_y <- c(out_y, y_val)
                            }
                          }

                          list(x = out_x, y = out_y)
                        },

                        .infer_scale_y = function() {
                           if (!is.null(self$gr_y)) {
                             if (self$style %in% c("triangle", "rectangle")) {
                               # Distance axis: continuous scale [0, yDistMax] in bp.
                               # yDistMax takes priority; fall back to maxDist.
                               dist_max <- self$yDistMax %||% self$maxDist
                               if (!is.null(dist_max))
                                 seq_scale_continuous(limits = c(0, dist_max))
                               else
                                 NULL  # data-range fallback handled in prep()
                             } else {
                               seq_scale_genomic(self$gr_y)
                             }
                           } else if (!is.null(self$groups)) {
                             seq_scale_discrete(levels = self$groups)
                           } else {
                             NULL
                           }
                         }
                       )
)

# Helper functions for SeqTile style support

#' @keywords internal
.validate_style_params <- function(style, gr_y, maxDist) {
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
    # Keep upper triangle only: y >= x (same as "diagonal").
    # Lower-triangle tiles produce negative rotated y-coordinates which
    # cause a symmetric filled rectangle instead of the expected triangular shape.
    y_start <- start(gr_y)[qh]
    x_start <- start(gr_x)[qh]
    return(y_start >= x_start)
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
