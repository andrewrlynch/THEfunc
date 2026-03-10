# ── Constructor wrappers ──────────────────────────────────────────────────────
# Rename R6 class bindings to hidden names so the public short-form constructors
# can share the original names.  The R6 class *string* (first arg to R6Class)
# is unchanged, so inherits(), class(), is(), etc. all continue to work.

.SeqTrackR6   <- SeqTrack
.SeqPlotR6    <- SeqPlot
.SeqPointR6   <- SeqPoint
.SeqTileR6    <- SeqTile
.SeqMatrixR6  <- SeqMatrix
.SeqSegmentR6 <- SeqSegment
.SeqRectR6    <- SeqRect
.SeqBarR6     <- SeqBar
.SeqLineR6    <- SeqLine
.SeqAreaR6    <- SeqArea
.SeqStringR6  <- SeqString
.SeqLinkR6    <- SeqLink
.SeqArchR6    <- SeqArch
.SeqIdeogramR6 <- SeqIdeogram
.SeqReconR6   <- SeqRecon
.SeqGeneR6    <- SeqGene

# Pin `inherit` fields to the stored R6ClassGenerators before the class names
# are overwritten with function wrappers.  R6 lazily evaluates `inherit =` at
# $new() time, so without this, classes that inherit a wrapped name (e.g.
# SeqArch$inherit → "SeqLink") would resolve to the function wrapper and error
# with "inherit must be a R6ClassGenerator".
#
# Affected chains:
#   SeqString      → SeqLink  (SeqLink will become a function)
#   SeqArch        → SeqLink
#   SeqRecon       → SeqArch  (SeqArch will become a function)
#   SeqHighlightClass → SeqLink  (not wrapped itself, but inherits SeqLink)
.SeqStringR6$inherit    <- .SeqLinkR6
.SeqArchR6$inherit      <- .SeqLinkR6
.SeqReconR6$inherit     <- .SeqArchR6
SeqHighlightClass$inherit <- .SeqLinkR6

#' Create a new SeqTrack
#' @inheritParams SeqTrack_R6
#' @export
SeqTrack <- function(...) .SeqTrackR6$new(...)

#' Create a new SeqPlot
#' @inheritParams SeqPlot_R6
#' @export
SeqPlot <- function(...) .SeqPlotR6$new(...)

#' @export
SeqPoint   <- function(...) .SeqPointR6$new(...)
#' @export
SeqTile    <- function(...) .SeqTileR6$new(...)
#' @export
SeqMatrix  <- function(...) .SeqMatrixR6$new(...)
#' @export
SeqSegment <- function(...) .SeqSegmentR6$new(...)
#' @export
SeqRect    <- function(...) .SeqRectR6$new(...)
#' @export
SeqBar     <- function(...) .SeqBarR6$new(...)
#' @export
SeqLine    <- function(...) .SeqLineR6$new(...)
#' @export
SeqArea    <- function(...) .SeqAreaR6$new(...)
#' @export
SeqString  <- function(...) .SeqStringR6$new(...)
#' @export
SeqLink    <- function(...) .SeqLinkR6$new(...)
#' @export
SeqArch    <- function(...) .SeqArchR6$new(...)
#' @export
SeqIdeogram <- function(...) .SeqIdeogramR6$new(...)
#' @export
SeqRecon   <- function(...) .SeqReconR6$new(...)
#' @export
SeqGene    <- function(...) .SeqGeneR6$new(...)

# ── Operator API ──────────────────────────────────────────────────────────────

#' Stack a SeqTrack into a SeqPlot
#'
#' The `%|%` operator adds a `SeqTrack` to a `SeqPlot`, making it the new
#' "current" track (subsequent `%+%` calls add elements to it).
#'
#' @param e1 A `SeqPlot` object.
#' @param e2 A `SeqTrack` object.
#' @return The `SeqPlot` (invisibly modified in place).
#' @export
`%|%` <- function(e1, e2) {
  if (inherits(e1, "SeqPlot")) {
    if (!inherits(e2, "SeqTrack"))
      stop("%|% expects a SeqTrack on the right-hand side (got '", class(e2)[1], "')")
    e1$addTrack(e2)
    return(e1)
  }
  stop("No method for '%|%' for object of class '", class(e1)[1], "'")
}

#' Add a SeqElement to a SeqTrack or SeqPlot
#'
#' When applied to a `SeqPlot`, `%+%` appends the element to the last (current)
#' track.  When applied to a `SeqTrack` directly, it appends the element to
#' that track.
#'
#' @param e1 A `SeqPlot` or `SeqTrack`.
#' @param e2 A `SeqElement` or `SeqLink`.
#' @return `e1`, invisibly modified in place.
#' @export
`%+%` <- function(e1, e2) {
  if (inherits(e1, "SeqPlot")) {
    n <- length(e1$tracks)
    if (n == 0) stop("No tracks in SeqPlot yet — start with SeqPlot() %|% SeqTrack()")
    if (inherits(e2, c("SeqElement", "SeqLink"))) {
      e1$tracks[[n]]$addElement(e2)
      return(e1)
    }
    if (inherits(e2, "SeqTrack")) {
      e1$addTrack(e2)
      return(e1)
    }
    stop("%+% cannot add object of class '", class(e2)[1], "' to a SeqPlot")
  }
  if (inherits(e1, "SeqTrack")) {
    if (!inherits(e2, c("SeqElement", "SeqLink")))
      stop("%+% cannot add object of class '", class(e2)[1], "' to a SeqTrack")
    e1$addElement(e2)
    return(e1)
  }
  stop("No method for '%+%' for object of class '", class(e1)[1], "'")
}
