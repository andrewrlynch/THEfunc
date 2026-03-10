# SeqLine ----
#' SeqLine R6 Class
#'
#' @description
#' R6 class for drawing line plots from genomic intervals in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' Each genomic interval is represented as a point at its midpoint by default,
#' with y-values taken from a metadata column or set to a constant. The points
#' are connected into a continuous line. Step lines can also be drawn by setting
#' the `type` aesthetic to `"s"` or `"step"`.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100, 200), width = 50),
#'   score = c(2, 5, 3)
#' )
#' line <- SeqLine(gr, yCol = "score")
#' line$prep(layout_track = some_layout, track_windows = some_windows)
#' line$draw()
#'
#' @export
SeqLine <- R6::R6Class("SeqLine",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic intervals.
                         gr = NULL,

                         #' @field y Numeric vector of y-values (from `yCol` or constant).
                         y = NULL,

                         #' @field yCol Optional column name in `gr` used for y-values.
                         yCol = NULL,

                         #' @field coordOriginal A `GRanges` object storing unmodified input coordinates.
                         coordOriginal = NULL,

                         #' @field coordCanvas List of matrices storing transformed line coordinates
                         #'   in canvas space for each genomic window.
                         coordCanvas = NULL,

                         #' @field aesthetics List of current aesthetics merged with defaults.
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default aesthetics for lines:
                         #'   \code{type = "n"}, \code{size = 0.1}, \code{color = "#1C1B1A"}.
                         defaultAesthetics = list(
                           type = "n",
                           size = 0.1,
                           color = "#1C1B1A"
                         ),

                         #' @description
                         #' Create a new `SeqLine` object.
                         #'
                         #' @param gr A `GRanges` object containing genomic intervals.
                         #' @param yCol Optional column name in `gr` for y-values.
                         #' @param aesthetics Optional list of aesthetic overrides.
                         #' @return A new `SeqLine` object.
                         initialize = function(gr = NULL, x = NULL, yCol = NULL, aesthetics = list(),
                                               aes = NULL, scale = NULL) {
                           if (!is.null(x)) {
                             stopifnot(inherits(x, "GRanges"))
                             gr_use <- x
                           } else {
                             stopifnot(inherits(gr, "GRanges"))
                             gr_use <- gr
                           }
                           self$gr <- gr_use
                           self$coordOriginal <- gr_use
                           self$yCol <- yCol

                           if (!is.null(yCol) && yCol %in% names(mcols(gr_use))) {
                             self$y <- as.numeric(mcols(gr_use)[[yCol]])
                           } else {
                             self$y <- rep(0.5, length(gr_use))
                           }

                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)

                           if (!is.null(aes)) {
                             resolved <- .resolve_aes(
                               data_mcols    = as.data.frame(S4Vectors::mcols(gr)),
                               aes_obj       = aes, scale_obj = scale,
                               n             = length(gr), default_color = "#1C1B1A"
                             )
                             for (nm in names(resolved)) self$aesthetics[[nm]] <- resolved[[nm]]
                           }
                         },

                         #' @description
                         #' Prepare line coordinates in canvas space for each genomic window.
                         #'
                         #' @param layout_track A list of panel layout metadata for a track.
                         #' @param track_windows A `GRanges` object defining genomic windows.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))

                           ov <- findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- S4Vectors::queryHits(ov)
                           sh <- S4Vectors::subjectHits(ov)

                           x <- (start(self$gr)[qh] + end(self$gr)[qh]) / 2
                           y <- self$y[qh]

                           if (self$aesthetics$type %in% c("s", "step")) {
                             x <- start(self$gr)[qh]
                             x <- rep(x, each = 2)[-1]
                             y <- rep(y, each = 2)[-length(y) * 2]
                           }

                           for (w in unique(sh)) {
                             mask <- sh == w
                             if (sum(mask) == 0) next

                             x_sub <- x[mask]
                             y_sub <- y[mask]
                             panel_meta <- layout_track[[w]]

                             u <- (x_sub - panel_meta$xscale[1]) / diff(panel_meta$xscale)
                             v <- (y_sub - panel_meta$yscale[1]) / diff(panel_meta$yscale)
                             u <- pmax(pmin(u, 1), 0)
                             v <- pmax(pmin(v, 1), 0)

                             x_canvas <- panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0)
                             y_canvas <- panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)

                             self$coordCanvas[[w]] <- cbind(x_canvas, y_canvas)
                           }

                           invisible()
                         },

                         #' @description
                         #' Draw lines onto the plotting canvas.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()
                           for (w in seq_along(self$coordCanvas)) {
                             coords <- self$coordCanvas[[w]]
                             if (is.null(coords)) next
                             grid.lines(
                               x = unit(coords[, 1], "npc"),
                               y = unit(coords[, 2], "npc"),
                               gp = gpar(col = self$aesthetics$color, cex = self$aesthetics$size)
                             )
                           }
                         }
                       )
)