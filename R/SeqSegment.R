# SeqSegment
#' SeqSegment R6 Class
#'
#' @description
#' R6 class for plotting line segments in the SeqPlot R6 framework.
#' Segments represent features with genomic start–end ranges and y-values
#' from metadata or defaults.
#'
#' @details
#' The `SeqSegment` class inherits from [SeqElement] and draws horizontal or
#' vertical line segments. Y-values can be supplied from a single column
#' (`yCol`) or from separate start and end columns (`y0Col`, `y1Col`).
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(
#'   "chr1",
#'   IRanges(c(1, 100), width = 50),
#'   score = c(0.2, 0.8)
#' )
#' seg <- SeqSegment(gr, yCol = "score")
#' seg$prep(layout_track = layout_info[[1]], track_windows = global_windows)
#' seg$draw()
#'
#' @export
SeqSegment <- R6::R6Class(
  "SeqSegment",
  inherit = SeqElement,
  public = list(

    #' @field gr A `GRanges` object containing genomic intervals.
    gr = NULL,

    #' @field y0 Numeric vector of lower y-values.
    y0 = NULL,

    #' @field y1 Numeric vector of upper y-values.
    y1 = NULL,

    #' @field yCol Optional metadata column for both y0 and y1 values.
    yCol = NULL,

    #' @field y0Col Optional metadata column for lower y-values.
    y0Col = NULL,

    #' @field y1Col Optional metadata column for upper y-values.
    y1Col = NULL,

    #' @field coordOriginal A `GRanges` object storing original input coordinates.
    coordOriginal = NULL,

    #' @field coordCanvas A list of transformed segment coordinates per window.
    coordCanvas = NULL,

    #' @field aesthetics List of visual aesthetics.
    aesthetics = NULL,

    #' @field defaultAesthetics Default values for aesthetics.
    defaultAesthetics = list(
      lwd = 1.5,
      col = "#1C1B1A"
    ),

    #' @description
    #' Create a new `SeqSegment` object.
    #'
    #' @param gr A `GRanges` object containing genomic intervals.
    #' @param yCol Column name in `gr` for both y0 and y1 values (optional).
    #' @param y0Col Column name in `gr` for lower y-values (optional).
    #' @param y1Col Column name in `gr` for upper y-values (optional).
    #' @param aesthetics Optional list of aesthetic overrides.
    initialize = function(gr = NULL, x = NULL, yCol = NULL, y0Col = NULL, y1Col = NULL,
                          aesthetics = list(), aes = NULL, scale = NULL) {
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
      self$y0Col <- y0Col
      self$y1Col <- y1Col

      # Parse y-values
      if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
        self$y0 <- self$y1 <- as.numeric(mcols(gr)[[yCol]])
      } else {
        if (!is.null(y0Col) && y0Col %in% names(mcols(gr))) {
          self$y0 <- as.numeric(mcols(gr)[[y0Col]])
        } else {
          self$y0 <- rep(0.5, length(gr))
        }
        if (!is.null(y1Col) && y1Col %in% names(mcols(gr))) {
          self$y1 <- as.numeric(mcols(gr)[[y1Col]])
        } else {
          self$y1 <- self$y0
        }
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
    #' Prepare segment coordinates by mapping genomic intervals into canvas space.
    #'
    #' @param layout_track Panel layout metadata (from grid layout).
    #' @param track_windows Genomic windows (`GRanges`) for this track.
    prep = function(layout_track, track_windows) {
      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(self$gr, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      x0 <- start(self$gr)[qh]
      x1 <- end(self$gr)[qh]
      y0 <- self$y0[qh]
      y1 <- self$y1[qh]

      for (w in unique(sh)) {
        mask <- sh == w
        if (sum(mask) == 0) next
        panel_meta <- layout_track[[w]]

        u0 <- (x0[mask] - panel_meta$xscale[1]) / diff(panel_meta$xscale)
        u1 <- (x1[mask] - panel_meta$xscale[1]) / diff(panel_meta$xscale)
        v0 <- (y0[mask] - panel_meta$yscale[1]) / diff(panel_meta$yscale)
        v1 <- (y1[mask] - panel_meta$yscale[1]) / diff(panel_meta$yscale)

        u0 <- pmax(pmin(u0, 1), 0)
        u1 <- pmax(pmin(u1, 1), 0)
        v0 <- pmax(pmin(v0, 1), 0)
        v1 <- pmax(pmin(v1, 1), 0)

        x0_canvas <- panel_meta$inner$x0 + u0 * (panel_meta$inner$x1 - panel_meta$inner$x0)
        x1_canvas <- panel_meta$inner$x0 + u1 * (panel_meta$inner$x1 - panel_meta$inner$x0)
        y0_canvas <- panel_meta$inner$y0 + v0 * (panel_meta$inner$y1 - panel_meta$inner$y0)
        y1_canvas <- panel_meta$inner$y0 + v1 * (panel_meta$inner$y1 - panel_meta$inner$y0)

        self$coordCanvas[[w]] <- list(
          x0 = x0_canvas,
          x1 = x1_canvas,
          y0 = y0_canvas,
          y1 = y1_canvas
        )
      }
      invisible()
    },

    #' @description
    #' Draw line segments on the plotting canvas.
    draw = function() {
      if (is.null(self$coordCanvas)) return()
      for (coords in self$coordCanvas) {
        if (is.null(coords)) next
        grid::grid.segments(
          x0 = grid::unit(coords$x0, "npc"),
          x1 = grid::unit(coords$x1, "npc"),
          y0 = grid::unit(coords$y0, "npc"),
          y1 = grid::unit(coords$y1, "npc"),
          gp = grid::gpar(
            fill = self$aesthetics$fill,
            col = self$aesthetics$col,
            lwd = self$aesthetics$lwd,
            lineend = "butt"
          )
        )
      }
    }
  )
)