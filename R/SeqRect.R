# SeqRect ----
#' SeqRect R6 Class
#'
#' @description
#' R6 class for drawing rectangular genomic features (boxes) in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' `SeqRect` is useful for visualizing genomic intervals as filled rectangles,
#' for example in bar plots, ideograms, or feature blocks. By default, each
#' rectangle is drawn centered on a y-value (`yCol` or a constant) with a
#' configurable relative height (`width` aesthetic).
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100), width = 50),
#'   score = c(0.2, 0.8)
#' )
#' rects <- SeqRect(gr, yCol = "score")
#' rects$prep(layout_track = some_layout, track_windows = some_windows)
#' rects$draw()
#'
#' @export
SeqRect <- R6::R6Class("SeqRect",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic intervals.
                         gr = NULL,

                         #' @field y0 Numeric vector of lower rectangle y-values (computed in prep).
                         y0 = NULL,

                         #' @field y1 Numeric vector of upper rectangle y-values (computed in prep).
                         y1 = NULL,

                         #' @field y Numeric vector of rectangle center y-values.
                         y = NULL,

                         #' @field yCol Optional column name in `gr` used for y-values.
                         yCol = NULL,

                         #' @field coordCanvas A list of matrices containing transformed rectangle
                         #'   coordinates for each window (x0, x1, y0, y1).
                         coordCanvas = NULL,

                         #' @field aesthetics List of current aesthetics merged with defaults.
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default aesthetics: \code{fill = "grey80"},
                         #'   \code{col = "#1C1B1A"}, \code{lwd = 0.5}, \code{width = 0.1}.
                         defaultAesthetics = list(
                           fill = "grey80",
                           col = "#1C1B1A",
                           lwd = 0.5,
                           width = 0.1
                         ),

                         #' @description
                         #' Create a new `SeqRect` object.
                         #'
                         #' @param gr A `GRanges` object containing genomic intervals.
                         #' @param yCol Optional column name in `gr` used for rectangle y-centers.
                         #' @param aesthetics Optional list of aesthetic overrides.
                         #' @return A new `SeqRect` object.
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
                           self$yCol <- yCol
                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)

                           y_center <- if (!is.null(yCol) && yCol %in% names(mcols(gr_use))) {
                             as.numeric(mcols(gr_use)[[yCol]])
                           } else {
                             rep(0.5, length(gr_use))
                           }

                           self$y <- y_center

                           if (!is.null(aes)) {
                             resolved <- .resolve_aes(
                               data_mcols    = as.data.frame(S4Vectors::mcols(gr_use)),
                               aes_obj       = aes, scale_obj = scale,
                               n             = length(gr_use), default_color = "#1C1B1A"
                             )
                             for (nm in names(resolved)) self$aesthetics[[nm]] <- resolved[[nm]]
                           }
                         },

                         #' @description
                         #' Prepare rectangle coordinates by mapping genomic intervals into
                         #' panel-relative canvas space.
                         #'
                         #' @param layout_track A list of panel layout metadata.
                         #' @param track_windows A `GRanges` object of genomic windows.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))
                           ov <- findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- queryHits(ov)
                           sh <- subjectHits(ov)
                           x0 <- start(self$gr)[qh]
                           x1 <- end(self$gr)[qh]
                           y_center <- self$y[qh]

                           for (w in unique(sh)) {
                             mask <- sh == w
                             p <- layout_track[[w]]

                             # Subset to current panel's hits
                             x0_sub <- x0[mask]
                             x1_sub <- x1[mask]
                             y_sub <- y_center[mask]

                             # Clip horizontally
                             clip <- clipToXscale(x0_sub, x1_sub, p$xscale)

                             if (length(clip$x0) == 0) {
                               self$coordCanvas[[w]] <- NULL
                               next
                             }

                             x0_sub <- clip$x0
                             x1_sub <- clip$x1
                             y_sub  <- y_sub[clip$mask]

                             # Sanity check
                             if (length(y_sub) == 0) {
                               self$coordCanvas[[w]] <- matrix(
                                 numeric(0),
                                 ncol = 4,
                                 dimnames = list(NULL, c("x0", "x1", "y0", "y1"))
                               )
                               next
                             }

                             # Transform X to panel-relative [0–1]
                             u0 <- (x0_sub - p$xscale[1]) / diff(p$xscale)
                             u1 <- (x1_sub - p$xscale[1]) / diff(p$xscale)

                             # Vertical position and height
                             v_center <- (y_sub - p$yscale[1]) / diff(p$yscale)
                             v_center <- pmax(pmin(v_center, 1), 0)
                             v_half <- self$aesthetics$width / 2
                             v0 <- v_center - v_half
                             v1 <- v_center + v_half

                             # Convert to canvas space
                             xleft <- p$inner$x0 + u0 * (p$inner$x1 - p$inner$x0)
                             xright <- p$inner$x0 + u1 * (p$inner$x1 - p$inner$x0)
                             ybottom <- p$inner$y0 + v0 * (p$inner$y1 - p$inner$y0)
                             ytop <- p$inner$y0 + v1 * (p$inner$y1 - p$inner$y0)

                             # Store as matrix for vectorized drawing
                             self$coordCanvas[[w]] <- cbind(
                               x0 = xleft,
                               x1 = xright,
                               y0 = ybottom,
                               y1 = ytop
                             )
                           }
                         },

                         #' @description
                         #' Draw rectangles on the plotting canvas.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()

                           for (coords in self$coordCanvas) {
                             if (!is.matrix(coords) || nrow(coords) == 0 || is.null(colnames(coords))) next

                             grid.rect(
                               x = unit((coords[, "x0"] + coords[, "x1"]) / 2, "npc"),
                               y = unit((coords[, "y0"] + coords[, "y1"]) / 2, "npc"),
                               width = unit(abs(coords[, "x1"] - coords[, "x0"]), "npc"),
                               height = unit(abs(coords[, "y1"] - coords[, "y0"]), "npc"),
                               gp = gpar(
                                 fill = self$aesthetics$fill,
                                 col = self$aesthetics$col,
                                 lwd = self$aesthetics$lwd
                               )
                             )
                           }
                         }
                       )
)