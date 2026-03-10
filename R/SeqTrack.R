# SeqTrack ----
#' SeqTrack R6 Class
#'
#' @description
#' An R6 class representing a single track in the SeqPlot framework.
#' A track contains one or more sequence elements (`SeqElement` objects),
#' associated genomic windows (`GRanges`), and track-specific aesthetics
#' such as axis visibility and titles.
#'
#' @details
#' Tracks are the core containers in a `SeqPlot` object. Each track defines
#' its own genomic windows and holds a list of plotting elements (such as
#' `SeqPoint`, `SeqBar`, `SeqLine`, or `SeqLink`). Track-level aesthetics
#' allow customization of axes, labels, and y-axis ranges.
#'
#' @examples
#' library(GenomicRanges)
#' win <- GRanges("chr1", IRanges(c(1, 1001), width = 500))
#' track <- SeqTrack(windows = win)
#' track$addElement(SeqRect(win))
#' length(track$elements)
#'
#' @export
SeqTrack <- R6::R6Class("SeqTrack",
                        public = list(

                          #' @field elements A list of `SeqElement` objects contained in this track.
                          elements = NULL,

                          #' @field windows A `GRanges` object defining the genomic windows
                          #'   covered by this track.
                          windows = NULL,

                          #' @field scale_x Position scale for the x-axis (`SeqPositionScale`
                          #'   or `NULL`). If `NULL`, auto-built from `windows`.
                          scale_x = NULL,

                          #' @field scale_y Position scale for the y-axis (`SeqPositionScale`
                          #'   or `NULL`). If `NULL`, auto-inferred from element data.
                          scale_y = NULL,

                          #' @field y_windows Optional `GRanges` defining genomic y-axis windows
                          #'   when the track uses a genomic y-axis (for 2D plots).
                          y_windows = NULL,
                          #' @field uses_genomic_y Logical flag indicating whether elements in
                          #'   this track use genomic y-axis (`gr_y` provided).
                          uses_genomic_y = FALSE,

                          #' @field aesthetics A named list of aesthetics controlling axis
                          #'   rendering and labeling:
                          #'   \describe{
                          #'     \item{xAxisTitle}{Logical; whether to show the x-axis title.}
                          #'     \item{xAxisTitleText}{Optional custom text for the x-axis title.}
                          #'     \item{yAxisTitle}{Logical; whether to show the y-axis title.}
                          #'     \item{yAxisTitleText}{Optional custom text for the y-axis title.}
                          #'     \item{yAxisLimits}{Optional numeric vector of length 2 setting the y-axis limits.}
                          #'     \item{yAxisTitleRotation}{Optional numeric vector of length 1 setting the y-axis title rotation. Default is 90.}
                          #'   }
                          aesthetics = list(
                            xAxisTitle = TRUE,
                            xAxisTitleText = TRUE,
                            yAxisTitle = TRUE,
                            yAxisTitleText = NULL,
                            yAxisLimits = NULL,
                            yAxisTitleRotation = NULL,
                            yAxisGenomicLabels = FALSE
                          ),

                          #' @description
                          #' Create a new `SeqTrack` object.
                          #'
                          #' @param elements Optional list of `SeqElement` objects to add
                          #'   when constructing the track.
                          #' @param windows Optional `GRanges` object defining the genomic
                          #'   windows for this track. Sugar for
                          #'   `scale_x = seq_scale_genomic(windows)`.
                          #' @param scale_x Optional position scale for the x-axis
                          #'   (`SeqPositionScale`).
                          #' @param scale_y Optional position scale for the y-axis
                          #'   (`SeqPositionScale`).
                          #' @param aesthetics Optional named list of track aesthetics,
                          #'   overriding defaults.
                          #' @return A new `SeqTrack` object.
                          initialize = function(elements = list(),
                                                windows = NULL,
                                                y_windows = NULL,
                                                scale_x = NULL,
                                                scale_y = NULL,
                                                aesthetics = list(
                                                  xAxisTitle = TRUE,
                                                  xAxisTitleText = TRUE,
                                                  yAxisTitle = TRUE,
                                                  yAxisTitleText = NULL,
                                                  yAxisLimits = NULL,
                                                  yAxisTitleRotation = NULL,
                                                  yAxisGenomicLabels = FALSE
                                                )) {
                            self$elements <- elements
                            self$aesthetics <- aesthetics

                            # windows = sugar for scale_x = seq_scale_genomic(windows)
                            if (!is.null(windows) && is.null(scale_x)) {
                              self$scale_x <- seq_scale_genomic(windows)
                            } else if (!is.null(scale_x)) {
                              self$scale_x <- scale_x
                            }

                            # Always store windows for backward compat
                            if (!is.null(windows)) {
                              self$windows <- windows
                            } else if (inherits(scale_x, "SeqScaleGenomic")) {
                              self$windows <- scale_x$windows
                            }

                            # y_windows sugar for genomic y-axis
                            if (!is.null(y_windows) && is.null(scale_y)) {
                              self$y_windows <- y_windows
                              self$scale_y <- seq_scale_genomic(y_windows)
                            } else if (!is.null(scale_y)) {
                              self$scale_y <- scale_y
                              if (inherits(scale_y, "SeqScaleGenomic")) self$y_windows <- scale_y$windows
                            }

                            # track-level flag: whether elements in this track use genomic y
                            self$uses_genomic_y <- any(vapply(elements, function(e) !is.null(e$gr_y), logical(1)))
                            self$elements <- elements
                            self$aesthetics <- aesthetics
                          },

                          #' @description
                          #' Add a `SeqElement` object to this track.
                          #'
                          #' @param feature A `SeqElement` object to be appended to the track.
                          #' @return Updates the `elements` field with the new feature.
                          addElement = function(feature) {
                            if (!inherits(feature, "SeqElement")) stop("feature must be a SeqElement")

                            # Only SeqTile may provide genomic ranges for both axes
                            if (!is.null(feature$gr_y) && !inherits(feature, "SeqTile")) {
                              stop("Only SeqTile may have GRanges on both x and y axes.")
                            }

                            # Enforce consistency: all elements must agree on whether
                            # they use genomic y-axis (i.e., provide `gr_y`) or not.
                            if (length(self$elements) > 0) {
                              existing_genomic_y <- self$uses_genomic_y
                              new_genomic_y <- !is.null(feature$gr_y)
                              if (existing_genomic_y != new_genomic_y) {
                                stop("All elements in a SeqTrack must either all use genomic y-axis (provide `y=` GRanges) or none do.")
                              }
                            }

                            self$elements <- append(self$elements, list(feature))
                            self$uses_genomic_y <- any(vapply(self$elements, function(e) !is.null(e$gr_y), logical(1)))
                          }
                        )
)



# ---- SeqBlankTrack
#' Create a blank (axisless) track
#'
#' Produces a SeqTrack with no elements and with aesthetics set so that the
#' track draws as a pure blank spacer (no axes, ticks, labels, titles, borders).
#'
#' @param windows Optional GRanges for this track (defaults to seqPlot$windows later).
#' @param height Optional numeric height for this track; if provided and the plot
#'   uses trackHeights, addBlankTrack() will append it.
#' @param aesthetics Optional list to override any of the blank defaults.
#' @return A SeqTrack
SeqBlankTrack <- function(windows = NULL,
                          height = NULL,
                          aesthetics = list(),
                          singleWindow = TRUE) {

  if (is.null(windows) && isTRUE(singleWindow)) {
    windows <- .SeqUnitWindow()
  }

  blank_aes <- list(
    trackBackground  = NA,
    trackBorder      = NA,
    windowBackground = NA,
    windowBorder     = NA,

    xAxisLine   = FALSE,
    yAxisLine   = FALSE,
    xAxisTicks  = FALSE,
    yAxisTicks  = FALSE,
    xAxisLabels = FALSE,
    yAxisLabels = FALSE,
    xAxisTitle  = FALSE,
    yAxisTitle  = FALSE,

    xAxisBreakLines = FALSE,
    yAxisBreakLines = FALSE
  )

  blank_aes <- modifyList(blank_aes, aesthetics)

  trk <- SeqTrack(
    elements = list(),
    windows = windows,
    aesthetics = blank_aes
  )

  attr(trk, "blank_track_height") <- height
  trk
}