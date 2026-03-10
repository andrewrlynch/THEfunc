# SeqHighlight ----
#' SeqHighlight R6 Class
#'
#' @description
#' Draw a filled polygon "connector" that highlights the same genomic interval
#' between two tracks (e.g. an ideogram track and a zoomed view track).
#'
#' @details
#' `SeqHighlight` projects each interval in `gr` onto track `t0` and `t1`,
#' then draws a quadrilateral between the two projected spans.
#'
#' @param gr A `GRanges` of intervals to highlight.
#' @param t0 Track index for the first projection.
#' @param t1 Track index for the second projection.
#' @param side Which edge of each track to attach to: "auto", "top", "bottom".
#'   If "auto", uses bottom edge of the upper track and top edge of the lower track.
#' @param fill Fill color.
#' @param alpha Fill alpha in [0,1].
#' @param col Border color (NA for no border).
#' @param lwd Border width.
#'
SeqHighlightClass <- R6::R6Class(
  "SeqHighlight",
  inherit = SeqLink,
  public = list(

    side = NULL,
    aesthetics = NULL,
    coordGrid = NULL,

    initialize = function(gr, t0 = 1, t1 = 2,
                          side = c("auto", "top", "bottom"),
                          fill = "black", alpha = 0.15,
                          col = NA, lwd = 0.5) {

      side <- match.arg(side)

      if (!inherits(gr, "GRanges")) stop("gr must be a GRanges object.")
      self$gr1 <- gr
      self$gr2 <- gr  # unused, but keeps SeqLink contract simple

      self$t0 <- rep.int(t0, length(gr))
      self$t1 <- rep.int(t1, length(gr))

      self$side <- side
      self$aesthetics <- list(fill = fill, alpha = alpha, col = col, lwd = lwd)
    },

    prep = function(layout_all_tracks, track_windows_list, arc_track_idx = NULL) {

      # helper: project GRanges intervals onto a given track's panel(s) in root-npc
      project_to_track <- function(gr, track_idx, attach_edge = c("top", "bottom")) {
        attach_edge <- match.arg(attach_edge)

        layout_track  <- layout_all_tracks[[track_idx]]
        track_windows <- track_windows_list[[track_idx]]

        ov <- GenomicRanges::findOverlaps(gr, track_windows)
        if (length(ov) == 0) return(NULL)

        qh <- S4Vectors::queryHits(ov)
        sh <- S4Vectors::subjectHits(ov)

        out <- vector("list", length(qh))
        j <- 0L

        for (k in seq_along(qh)) {
          i <- qh[k]            # which gr row
          w <- sh[k]            # which window/panel
          p <- layout_track[[w]]

          x0 <- GenomicRanges::start(gr)[i]
          x1 <- GenomicRanges::end(gr)[i]

          # clip to the panel xscale (same strategy as SeqRect)
          clip <- clipToXscale(x0, x1, p$xscale)
          if (length(clip$x0) == 0) next

          x0c <- clip$x0
          x1c <- clip$x1

          u0 <- (x0c - p$xscale[1]) / diff(p$xscale)
          u1 <- (x1c - p$xscale[1]) / diff(p$xscale)

          xleft  <- p$full$x0 + u0 * (p$full$x1 - p$full$x0)
          xright <- p$full$x0 + u1 * (p$full$x1 - p$full$x0)

          yedge <- if (attach_edge == "top") p$full$y1 else p$full$y0

          j <- j + 1L
          out[[j]] <- data.frame(
            idx = i,
            track = track_idx,
            window = w,
            x0 = xleft,
            x1 = xright,
            y = yedge,
            stringsAsFactors = FALSE
          )
        }

        out <- out[seq_len(j)]
        if (length(out) == 0) return(NULL)
        do.call(rbind, out)
      }

      # decide attachment edges
      t0 <- self$t0[1]
      t1 <- self$t1[1]

      if (self$side == "auto") {
        if (t0 < t1) {
          edge0 <- "bottom"  # upper track attaches from its bottom
          edge1 <- "top"     # lower track attaches from its top
        } else if (t0 > t1) {
          edge0 <- "top"
          edge1 <- "bottom"
        } else {
          # same track: default to top
          edge0 <- "top"
          edge1 <- "top"
        }
      } else {
        edge0 <- self$side
        edge1 <- self$side
      }

      a <- project_to_track(self$gr1, t0, edge0)
      b <- project_to_track(self$gr1, t1, edge1)

      if (is.null(a) || is.null(b)) {
        self$coordGrid <- NULL
        return(invisible())
      }

      # pair by original GRanges row index
      # (assumes each interval typically hits one window per track; if multiple hits occur,
      #  this will draw one polygon per matching idx combination)
      merged <- merge(a, b, by = "idx", suffixes = c("0", "1"))
      if (nrow(merged) == 0) {
        self$coordGrid <- NULL
        return(invisible())
      }

      # store polygon corners
      self$coordGrid <- transform(
        merged,
        x00 = x00, x10 = x10, y0 = y0,
        x01 = x01, x11 = x11, y1 = y1
      )
      invisible()
    },

    draw = function() {
      df <- self$coordGrid
      if (is.null(df) || nrow(df) == 0) return()

      fill <- self$aesthetics$fill %||% "black"
      alpha <- self$aesthetics$alpha %||% 0.15
      col <- self$aesthetics$col %||% NA
      lwd <- self$aesthetics$lwd %||% 0.5

      fill2 <- grDevices::adjustcolor(fill, alpha.f = alpha)

      hlStemOffsetNPC = 0.01

      for (i in seq_len(nrow(df))) {
        xs <- c(df$x00[i], df$x00[i], df$x10[i], df$x10[i], df$x11[i], df$x11[i], df$x01[i], df$x01[i])
        ys <- c(df$y0[i] - hlStemOffsetNPC, df$y0[i], df$y0[i], df$y0[i] - hlStemOffsetNPC, df$y1[i] + hlStemOffsetNPC, df$y1[i],  df$y1[i], df$y1[i] + hlStemOffsetNPC)

        grid::grid.polygon(
          x = grid::unit(xs, "npc"),
          y = grid::unit(ys, "npc"),
          gp = grid::gpar(fill = fill2, col = col, lwd = lwd)
        )
      }
    }
  )
)

# Convenience constructor (matches your desired call style)
#' @export
SeqHighlight <- function(gr, t0 = 1, t1 = 2, ...) {
  SeqHighlightClass$new(gr = gr, t0 = t0, t1 = t1, ...)
}