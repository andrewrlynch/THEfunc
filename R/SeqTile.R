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
                         initialize = function(gr = NULL, x = NULL, y = NULL,
                                               groupCol = NULL, aesthetics = list(),
                                               aes = NULL, scale = NULL) {
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
                         },

                         #' @description
                         #' Prepare tile coordinates for plotting in canvas space.
                         #' @param layout_track A list of panel metadata for the current track.
                         #' @param track_windows A `GRanges` of genomic windows for the track.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))

                           ov <- GenomicRanges::findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- S4Vectors::queryHits(ov)
                           sh <- S4Vectors::subjectHits(ov)

                           x0     <- start(self$gr)[qh]
                           x1     <- end(self$gr)[qh]
                           colors <- as.character(mcols(self$gr)$color)[qh]

                           # If genomic y provided, compute y0/y1 from gr_y
                           y_idx <- NULL
                           y0_raw <- NULL
                           y1_raw <- NULL
                           if (!is.null(self$gr_y)) {
                             y0_raw <- start(self$gr_y)[qh]
                             y1_raw <- end(self$gr_y)[qh]
                           } else if (!is.null(self$y)) {
                             y_idx <- self$y[qh]
                           }

                           for (w in unique(sh)) {
                             mask <- sh == w
                             if (sum(mask) == 0) next

                             panel <- layout_track[[w]]

                             # x: genomic mapping
                             u0 <- (x0[mask] - panel$xscale[1]) / diff(panel$xscale)
                             u1 <- (x1[mask] - panel$xscale[1]) / diff(panel$xscale)
                             u0 <- pmax(pmin(u0, 1), 0)
                             u1 <- pmax(pmin(u1, 1), 0)

                             if (!is.null(y0_raw)) {
                               # Genomic y-axis
                               y0_m <- y0_raw[mask]
                               y1_m <- y1_raw[mask]

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
                                     col = colors[mask][y_in_sub]
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
                                   col = colors[mask]
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
                                 col = colors[mask]
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