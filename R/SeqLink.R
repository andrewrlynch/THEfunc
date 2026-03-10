# SeqLink ----
#' SeqLink R6 Class
#'
#' @description
#' R6 class representing a link between two genomic loci in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' A `SeqLink` connects two sets of genomic intervals (`gr1` and `gr2`),
#' optionally across different tracks. It stores track indices (`t0`, `t1`),
#' vertical anchor positions (`y0`, `y1`), and link orientations. Higher-level
#' link classes such as `SeqArch` extend this class to implement specific
#' visualization styles (arches, lines, bezier curves, etc.).
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 50), width = 1))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), width = 1))
#' link <- SeqLink(gr1, gr2, t0 = 1, t1 = 1, y0 = 0, y1 = 0, orientation = "+")
#' link$gr1
#'
#' @export
SeqLink <- R6::R6Class("SeqLink",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr1 A `GRanges` object containing the start loci of each link.
                         gr1 = NULL,

                         #' @field gr2 A `GRanges` object containing the end loci of each link.
                         gr2 = NULL,

                         #' @field t0 Integer vector giving the track indices of the start loci.
                         t0 = NULL,

                         #' @field t1 Integer vector giving the track indices of the end loci.
                         t1 = NULL,

                         #' @field y0 Numeric vector of vertical anchor positions for the start loci.
                         y0 = NULL,

                         #' @field y1 Numeric vector of vertical anchor positions for the end loci.
                         y1 = NULL,

                         #' @field orientation Character vector specifying orientation of each link
                         #'   (e.g. `"+"`, `"-"`, `"*"`).
                         orientation = NULL,

                         #' @description
                         #' Create a new `SeqLink` object.
                         #'
                         #' @param gr1 Deprecated start positions; use `x1` instead.
                         #' @param gr2 Deprecated end positions; use `x2` instead.
                         #' @param x1 Preferred: `GRanges` object giving start positions of links.
                         #' @param x2 Preferred: `GRanges` object giving end positions of links.
                         #' @param t0 Integer or vector giving track indices for start loci (default: `0`).
                         #' @param t1 Integer or vector giving track indices for end loci (default: `0`).
                         #' @param y0 Numeric or vector giving vertical start positions (default: `0`).
                         #' @param y1 Numeric or vector giving vertical end positions (default: `0`).
                         #' @param color Color for links (currently unused in base class).
                         #' @param orientation Character or vector specifying link orientation
                         #'   (default: `"*"`).
                         #' @return A new `SeqLink` object.
                         initialize = function(gr1 = NULL, gr2 = NULL, x1 = NULL, x2 = NULL, t0 = 0, t1 = 0,
                                               y0 = 0, y1 = 0,
                                               color = "black", orientation = "*") {
                           # Accept x1/x2 (preferred) or legacy gr1/gr2
                           if (!is.null(x1)) {
                             gr1_use <- x1
                           } else {
                             gr1_use <- gr1
                           }
                           if (!is.null(x2)) {
                             gr2_use <- x2
                           } else {
                             gr2_use <- gr2
                           }
                           stopifnot(length(gr1_use) == length(gr2_use))
                           self$gr1 <- gr1_use
                           self$gr2 <- gr2_use
                           self$t0 <- rep(t0, length(gr1_use))
                           self$t1 <- rep(t1, length(gr1_use))
                           self$y0 <- rep(y0, length(gr1_use))
                           self$y1 <- rep(y1, length(gr2_use))
                           self$orientation <- rep(orientation, length(gr1_use))
                         }
                       )
)