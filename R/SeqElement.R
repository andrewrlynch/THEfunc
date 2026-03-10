#' SeqElement R6 Class
#'
#' @description
#' Base R6 class representing a generic sequence element in the SeqPlot
#' framework. Stores the original genomic ranges and provides slots for
#' transformed coordinates in grid space.
#'
#' @details
#' The `SeqElement` class is the root R6 class for all genomic elements
#' in the plotting system. It is initialized with a `GRanges` object and
#' retains both the original genomic coordinates (`coordOriginal`) and
#' transformed grid coordinates (`coordGrid`).
#'
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:5, width = 1))
#' el <- SeqElement$new(gr)
#' el$gr
#'
#' @export
SeqElement <- R6::R6Class("SeqElement",
                          public = list(

                            #' @field gr A `GRanges` object containing the genomic coordinates
                            #'   (x-axis by default) and any associated metadata.
                            gr = NULL,

                            #' @field gr_y Optional `GRanges` when the element provides genomic
                            #'   coordinates for the y-axis (enables genomic y-axis / 2D features).
                            gr_y = NULL,

                            #' @field coordOriginal A `GRanges` object storing the unmodified input
                            #'   coordinates for the primary (x) axis.
                            coordOriginal = NULL,

                            #' @field coordGrid Placeholder for transformed coordinates after applying
                            #'   the SeqPlot grid layout.
                            coordGrid = NULL,

                            #' @description
                            #' Create a new `SeqElement` object.
                            #' Backwards-compatible: accepts either `gr=` or `x=` for the primary
                            #' genomic coordinate. An optional `y=` GRanges may be provided to
                            #' enable genomic y-axis behavior for compatible elements.
                            #' @param gr Deprecated. A `GRanges` object (old name for x).
                            #' @param x Preferred primary `GRanges` for x-axis positions.
                            #' @param y Optional `GRanges` for genomic y-axis positions.
                            #' @param yCol Optional character string naming a metadata column in `gr`
                            #'   to be used for numeric/categorical y-axis values.
                            #' @return A new `SeqElement` object.
                            initialize = function(gr = NULL, x = NULL, y = NULL, yCol = NULL) {
                              # Accept either x or gr for backwards compatibility
                              if (!is.null(x)) {
                                if (!inherits(x, "GRanges")) stop("x must be a GRanges object.")
                                self$gr <- x
                                self$coordOriginal <- x
                              } else if (!is.null(gr)) {
                                if (!inherits(gr, "GRanges")) stop("gr must be a GRanges object.")
                                self$gr <- gr
                                self$coordOriginal <- gr
                              } else {
                                stop("Either `x` or `gr` must be supplied and be a GRanges.")
                              }

                              if (!is.null(y)) {
                                if (!inherits(y, "GRanges")) stop("y must be a GRanges object when provided.")
                                self$gr_y <- y
                              }
                            },

                            #' @description
                            #' Infer the y-axis position scale from this element's data.
                            #' Subclasses override to return a `SeqPositionScale` when
                            #' the element dictates the y-axis type (e.g. SeqTile → discrete,
                            #' SeqMatrix → genomic).
                            #' @return A `SeqPositionScale` object, or `NULL`.
                            .infer_scale_y = function() NULL
                          )
)