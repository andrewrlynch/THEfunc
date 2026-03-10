# SeqMatrix ----
#' SeqMatrix R6 Class
#'
#' @description
#' An R6 class for plotting 2D genomic heatmaps (e.g. Hi-C).
#' Genomic ranges define both the x- and y-axes. Each pair of
#' bins is drawn as a rectangle, filled according to a value
#' matrix mapped through a color function.
#'
#' @export
SeqMatrix <- R6::R6Class("SeqMatrix",
                         inherit = SeqElement,
                         public = list(

                           #' @field gr_x A `GRanges` object for x-axis bins.
                           gr_x = NULL,

                           #' @field gr_y A `GRanges` object for y-axis bins.
                           gr_y = NULL,

                           #' @field value A numeric matrix of size length(gr_y) × length(gr_x).
                           value = NULL,

                           #' @field col_fun A function mapping numeric values to colors
                           #'   (e.g. from circlize::colorRamp2).
                           col_fun = NULL,

                           #' @field yCol y column
                           yCol = NULL,

                           #' @field y_vals y column
                           y_vals = NULL,

                           #' @field coordOriginal List storing original GRanges and values.
                           coordOriginal = NULL,

                           #' @field coordCanvas List of per-panel cell rectangles (x0, x1, y0, y1, color).
                           coordCanvas = NULL,

                           #' @field aesthetics Plot aesthetics (border, lwd).
                           aesthetics = NULL,

                           #' @field defaultAesthetics Default aesthetics for tiles.
                           defaultAesthetics = list(
                             border = NA,
                             lwd = 0.1
                           ),

                           #' @description
                           #' Create a new `SeqMatrix` object.
                           #' @param gr_x GRanges for x-axis bins.
                           #' @param gr_y GRanges for y-axis bins (if NULL, use gr_x).
                           #' @param value Numeric matrix (nrow = length(gr_y), ncol = length(gr_x)).
                           #' @param col_fun Function mapping numeric values to colors.
                           #' @param yCol y column
                           #' @param y_vals y values
                           #' @param aesthetics List of aesthetics (border, lwd).
                           initialize = function(gr_x, gr_y = NULL, value, col_fun,
                                                 yCol = NULL, y_vals = NULL, aesthetics = list()) {
                             stopifnot(inherits(gr_x, "GRanges"))
                             if (is.null(gr_y)) gr_y <- gr_x
                             stopifnot(inherits(gr_y, "GRanges"))

                             if (!is.matrix(value)) stop("`value` must be a matrix.")
                             if (nrow(value) != length(gr_y) || ncol(value) != length(gr_x)) {
                               stop("`value` must have dimensions length(gr_y) × length(gr_x).")
                             }
                             if (!is.function(col_fun)) stop("`col_fun` must be a function (e.g. circlize::colorRamp2).")

                             self$gr_x <- gr_x
                             self$gr_y <- gr_y
                             self$value <- value
                             self$col_fun <- col_fun
                             self$yCol <- yCol

                             if (!is.null(yCol) && yCol %in% names(mcols(gr_y))) {
                               self$y_vals <- as.numeric(mcols(gr_y)[[yCol]])
                             } else {
                               # default: y-axis genomic, same as x
                               self$y_vals <- c(min(start(gr_y)), max(end(gr_y)))
                             }

                             self$coordOriginal <- list(gr_x = gr_x, gr_y = gr_y, value = value)
                             self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)
                           },

                           #' @description
                           #' Prepare heatmap cell coordinates for plotting in canvas space.
                           #' @param layout_track List of panel metadata for the current track.
                           #' @param track_windows A `GRanges` of genomic windows for the track.
                           prep = function(layout_track, track_windows) {
                             self$coordCanvas <- vector("list", length(track_windows))

                             # Always use genomic y-scale if yCol is NULL
                             yscale <- if (!is.null(self$yCol) && self$yCol %in% names(mcols(self$gr_y))) {
                               range(self$y_vals, na.rm = TRUE)
                             } else {
                               c(min(start(self$gr_y)), max(end(self$gr_y)))
                             }

                             for (w in seq_along(track_windows)) {
                               # Note: yscale is now set by layoutGrid() via .infer_scale_y().
                               # Fall back to genomic range only if layoutGrid didn't set it.
                               if (is.null(layout_track[[w]]$y_scale_type) ||
                                   layout_track[[w]]$y_scale_type == "continuous") {
                                 layout_track[[w]]$yscale <- c(min(start(self$gr_y)), max(end(self$gr_y)))
                               }
                               panel_meta <- layout_track[[w]]

                               # coords
                               x0 <- start(self$gr_x)
                               x1 <- end(self$gr_x)
                               y0 <- start(self$gr_y)
                               y1 <- end(self$gr_y)

                               # normalize to [0,1] relative to panel_meta scales
                               u0 <- (x0 - panel_meta$xscale[1]) / diff(panel_meta$xscale)
                               u1 <- (x1 - panel_meta$xscale[1]) / diff(panel_meta$xscale)
                               v0 <- (y0 - panel_meta$yscale[1]) / diff(panel_meta$yscale)
                               v1 <- (y1 - panel_meta$yscale[1]) / diff(panel_meta$yscale)

                               # clip
                               u0 <- pmax(pmin(u0, 1), 0)
                               u1 <- pmax(pmin(u1, 1), 0)
                               v0 <- pmax(pmin(v0, 1), 0)
                               v1 <- pmax(pmin(v1, 1), 0)

                               # canvas coords
                               x0_canvas <- panel_meta$inner$x0 + u0 * (panel_meta$inner$x1 - panel_meta$inner$x0)
                               x1_canvas <- panel_meta$inner$x0 + u1 * (panel_meta$inner$x1 - panel_meta$inner$x0)
                               y0_canvas <- panel_meta$inner$y0 + v0 * (panel_meta$inner$y1 - panel_meta$inner$y0)
                               y1_canvas <- panel_meta$inner$y0 + v1 * (panel_meta$inner$y1 - panel_meta$inner$y0)

                               # expand value matrix into long form
                               df <- expand.grid(x = seq_along(x0), y = seq_along(y0))
                               df$x0 <- x0_canvas[df$x]
                               df$x1 <- x1_canvas[df$x]
                               df$y0 <- y0_canvas[df$y]
                               df$y1 <- y1_canvas[df$y]
                               df$val <- as.vector(self$value)
                               df$col <- self$col_fun(df$val)

                               self$coordCanvas[[w]] <- df

                               # print(self$coordCanvas[[w]])
                             }
                           },

                           #' @description
                           #' Draw the heatmap tiles.
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
                             if (!is.null(self$gr_y))
                               seq_scale_genomic(self$gr_y)
                             else
                               NULL
                           }
                         )
)