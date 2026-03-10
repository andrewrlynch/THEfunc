# SeqIdeogram ----
#' SeqIdeogram R6 Class
#'
#' @description
#' An R6 class for drawing chromosome ideograms from cytogenetic band
#' annotations. Each band is drawn as a rectangle, and centromeres
#' (`acen` regions) are rendered as paired red triangles.
#'
#' @export
SeqIdeogram <- R6::R6Class("SeqIdeogram",
                           inherit = SeqElement,
                           public = list(

                             #' @field cytobands A `GRanges` object of cytogenetic bands. Must
                             #'   contain a `gieStain` metadata column.
                             cytobands = NULL,

                             #' @field coordCanvas A list of data frames holding canvas
                             #'   coordinates for non-centromeric bands (`x0`, `x1`, `y0`, `y1`)
                             #'   and their fill colors.
                             coordCanvas = NULL,

                             #' @field centroPolys A list of polygons representing centromere
                             #'   regions, each stored as x/y coordinate vectors.
                             centroPolys = NULL,

                             #' @description
                             #' Create a new `SeqIdeogram` object.
                             #' @param cytobands A `GRanges` object of cytobands with a `gieStain`
                             #'   metadata column.
                             #' @return A new `SeqIdeogram` object.
                             #' @examples
                             #' cb <- GenomicRanges::GRanges(
                             #'   seqnames = Rle("chr1"),
                             #'   ranges   = IRanges::IRanges(c(1, 500001, 1000001), width = 500000),
                             #'   gieStain = c("gneg", "acen", "gpos100")
                             #' )
                             #' ideogram <- SeqIdeogram(cb)
                             initialize = function(cytobands) {
                               stopifnot(inherits(cytobands, "GRanges"))
                               self$cytobands <- cytobands
                             },

                             #' @description
                             #' Prepare ideogram band and centromere coordinates for each panel
                             #' window in the SeqPlot layout.
                             #' @param layout_track A list of panel metadata from
                             #'   `SeqPlot$layoutGrid()`.
                             #' @param track_windows A `GRanges` object of genomic windows for the
                             #'   track.
                             prep = function(layout_track, track_windows) {
                               n_panels <- length(track_windows)
                               self$coordCanvas <- vector("list", n_panels)
                               self$centroPolys <- vector("list", n_panels)

                               ov <- findOverlaps(self$cytobands, track_windows)
                               if (length(ov) == 0) return()

                               qh <- queryHits(ov); sh <- subjectHits(ov)

                               for (w in unique(sh)) {
                                 panel <- layout_track[[w]]
                                 bands <- self$cytobands[qh[sh == w]]

                                 u0 <- (start(bands) - panel$xscale[1]) / diff(panel$xscale)
                                 u1 <- (end(bands)   - panel$xscale[1]) / diff(panel$xscale)
                                 u0 <- pmax(pmin(u0, 1), 0); u1 <- pmax(pmin(u1, 1), 0)

                                 x0c <- panel$inner$x0 + u0 * (panel$inner$x1 - panel$inner$x0)
                                 x1c <- panel$inner$x0 + u1 * (panel$inner$x1 - panel$inner$x0)
                                 y0c <- panel$inner$y0
                                 y1c <- panel$inner$y1

                                 stain <- mcols(bands)$gieStain
                                 fillCols <- sapply(stain, function(s) {
                                   if (s == "gneg") {
                                     "#FFFFFF"
                                     }
                                   else if (startsWith(s, "gpos")) {
                                     pct <- 1-as.numeric(sub("gpos", "", s)) / 100
                                     grey(pct)
                                   } else if (s == "acen"){
                                     "#FF0000"
                                     }
                                   else if (s == "stalk"){
                                     "#7AC0CF"
                                   }
                                   else if (s == "gvar"){
                                     "#CCF5FF"
                                   } else {
                                     "#CCCCCC"
                                     }
                                 })

                                 non_acen <- stain != "acen"
                                 self$coordCanvas[[w]] <- data.frame(
                                   x0 = x0c[non_acen], x1 = x1c[non_acen],
                                   y0 = y0c,          y1 = y1c,
                                   fill = fillCols[non_acen],
                                   stringsAsFactors = FALSE
                                 )

                                 acen_idx <- which(stain == "acen")
                                 if (length(acen_idx) == 2) {
                                   x0a <- x0c[acen_idx[1]]; x1a <- x1c[acen_idx[1]]
                                   x0b <- x0c[acen_idx[2]]; x1b <- x1c[acen_idx[2]]
                                   ym  <- (y0c + y1c) / 2

                                   tri1_x <- c(x0a, x1a, x0a)
                                   tri1_y <- c(y0c, ym,  y1c)

                                   tri2_x <- c(x1b, x0b, x1b)
                                   tri2_y <- c(y0c, ym,  y1c)

                                   self$centroPolys[[w]] <- list(
                                     list(x = tri1_x, y = tri1_y),
                                     list(x = tri2_x, y = tri2_y)
                                   )
                                 }
                               }
                             },

                             #' @description
                             #' Draw the ideogram for all windows: non-centromeric bands as
                             #' rectangles, centromeres as paired red triangles.
                             draw = function() {
                               if (!is.null(self$coordCanvas)) {
                                 for (coords in self$coordCanvas) {
                                   if (nrow(coords) == 0) next
                                   grid.rect(
                                     x      = unit((coords$x0 + coords$x1) / 2, "npc"),
                                     y      = unit((coords$y0 + coords$y1) / 2, "npc"),
                                     width  = unit(coords$x1 - coords$x0, "npc"),
                                     height = unit(coords$y1 - coords$y0, "npc"),
                                     gp     = gpar(fill = coords$fill, col = "black", lwd = 0.1)
                                   )
                                 }
                               }
                               if (!is.null(self$centroPolys)) {
                                 for (polys in self$centroPolys) {
                                   if (is.null(polys)) next
                                   for (tri in polys) {
                                     grid.polygon(
                                       x = unit(tri$x, "npc"),
                                       y = unit(tri$y, "npc"),
                                       gp = gpar(fill = "#FF0000", col = "black", lwd = 0.1)
                                     )
                                   }
                                 }
                               }
                             }
                           )
)