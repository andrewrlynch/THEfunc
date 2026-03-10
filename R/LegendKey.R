# ---- LegendKey ---------------------------------------------------------------

LegendKey <- function(label = NULL,
                      color = "#1C1B1A",
                      shape = "-",
                      size = 1,
                      alpha = 1,
                      fill = NULL,
                      lty = 1,
                      ...) {
  key <- list(
    label = label,
    color = color,
    shape = shape,
    size  = size,
    alpha = alpha,
    fill  = fill,
    lty   = lty,
    extra = list(...)
  )
  class(key) <- "LegendKey"
  key
}

.SeqLegend_drawKey <- function(key, x0, x1, y, height = 0.04) {
  if (is.null(key) || !inherits(key, "LegendKey")) return(invisible())

  col   <- key$color
  alpha <- key$alpha
  lwd   <- key$size
  lty   <- key$lty

  if (identical(key$shape, "rect")) {
    fill <- key$fill
    if (is.null(fill)) fill <- col
    grid::grid.rect(
      x = grid::unit((x0 + x1) / 2, "npc"),
      y = grid::unit(y, "npc"),
      width  = grid::unit(max(0, x1 - x0), "npc"),
      height = grid::unit(height, "npc"),
      gp = grid::gpar(col = col, fill = fill, alpha = alpha, lwd = lwd)
    )
    return(invisible())
  }

  if (identical(key$shape, "-") || identical(key$shape, "line")) {
    grid::grid.lines(
      x = grid::unit(c(x0, x1), "npc"),
      y = grid::unit(c(y,  y),  "npc"),
      gp = grid::gpar(col = col, alpha = alpha, lwd = lwd, lty = lty)
    )
    return(invisible())
  }

  # point pch
  pch <- key$shape
  grid::grid.points(
    x = grid::unit((x0 + x1) / 2, "npc"),
    y = grid::unit(y, "npc"),
    pch = pch,
    gp = grid::gpar(col = col, alpha = alpha),
    size = grid::unit(key$size, "char")
  )

  invisible()
}

SeqLegendElement <- R6::R6Class(
  "SeqLegendElement",
  public = list(
    entries = NULL,
    ncol = NULL,
    nrow = NULL,
    colSpacing = 0,
    rowSpacing = 0,
    hjust = 0,
    vjust = 0.5,
    keyWidth = 0.04,
    keyGap = 0.006,
    rowHeight = 0.12,
    textCex = 0.9,
    textCol = "#1C1B1A",
    fontface = 1,
    coordCanvas = NULL,   # SINGLE layout: list(cells=...)

    initialize = function(...,
                          ncol = NULL, nrow = NULL,
                          colSpacing = 0, rowSpacing = 0,
                          hjust = 0, vjust = 0.5,
                          keyWidth = 0.04, keyGap = 0.006,
                          rowHeight = 0.12,
                          textCex = 0.9, textCol = "#1C1B1A",
                          fontface = 1) {

      items <- list(...)
      if (length(items) == 0) stop("SeqLegend needs at least one LegendKey().")

      nm <- names(items)
      for (i in seq_along(items)) {
        if (!inherits(items[[i]], "LegendKey")) stop("All SeqLegend entries must be LegendKey().")

        if (is.null(items[[i]]$label) || !nzchar(items[[i]]$label)) {
          if (!is.null(nm) && nzchar(nm[[i]])) {
            items[[i]]$label <- nm[[i]]
          } else {
            stop("Each LegendKey needs label=... or be named 'Label'=LegendKey(...).")
          }
        }
      }

      self$entries <- items
      self$ncol <- ncol
      self$nrow <- nrow
      self$colSpacing <- colSpacing
      self$rowSpacing <- rowSpacing
      self$hjust <- hjust
      self$vjust <- vjust
      self$keyWidth <- keyWidth
      self$keyGap <- keyGap
      self$rowHeight <- rowHeight
      self$textCex <- textCex
      self$textCol <- textCol
      self$fontface <- fontface
    },

    # We intentionally compute ONE layout (legend track has one window)
    prep = function(layout_track, track_windows) {
      n_items <- length(self$entries)
      if (n_items == 0) return(invisible())

      ncol <- self$ncol
      nrow <- self$nrow

      if (is.null(ncol) && is.null(nrow)) {
        nrow <- 1
        ncol <- n_items
      } else if (is.null(ncol)) {
        ncol <- ceiling(n_items / nrow)
      } else if (is.null(nrow)) {
        nrow <- ceiling(n_items / ncol)
      }

      ncol <- max(1L, as.integer(ncol))
      nrow <- max(1L, as.integer(nrow))

      pad_x <- 0.02
      pad_y <- 0.02

      total_col_gap <- self$colSpacing * max(0, ncol - 1)
      total_row_gap <- self$rowSpacing * max(0, nrow - 1)

      block_height <- min(1 - 2 * pad_y, nrow * self$rowHeight + total_row_gap)
      cell_h <- (block_height - total_row_gap) / nrow

      block_width <- 1 - 2 * pad_x
      cell_w <- (block_width - total_col_gap) / ncol

      x_left   <- pad_x + (block_width - (ncol * cell_w + total_col_gap)) * self$hjust
      y_bottom <- pad_y + ((1 - 2 * pad_y) - block_height) * (1 - self$vjust)

      cells <- vector("list", n_items)
      for (i in seq_len(n_items)) {
        r <- ((i - 1) %/% ncol) + 1
        c <- ((i - 1) %%  ncol) + 1

        x0 <- x_left + (c - 1) * (cell_w + self$colSpacing)
        x1 <- x0 + cell_w
        y1 <- y_bottom + block_height - (r - 1) * (cell_h + self$rowSpacing)
        y0 <- y1 - cell_h

        key <- self$entries[[i]]
        xk1 <- min(x1, x0 + self$keyWidth)

        cells[[i]] <- list(
          label  = key$label,
          key    = key,
          key_x0 = x0,
          key_x1 = xk1,
          y      = (y0 + y1) / 2,
          text_x = min(0.98, xk1 + self$keyGap),
          text_y = (y0 + y1) / 2
        )
      }

      self$coordCanvas <- list(cells = cells)
      invisible()
    },

    # Draw ONCE in the current viewport (no window loop!)
    draw = function() {
      if (is.null(self$coordCanvas) || is.null(self$coordCanvas$cells)) return(invisible())
      cells <- self$coordCanvas$cells
      if (length(cells) == 0) return(invisible())

      for (cell in cells) {
        .SeqLegend_drawKey(
          key = cell$key,
          x0  = cell$key_x0,
          x1  = cell$key_x1,
          y   = cell$y,
          height = min(0.06, self$rowHeight * 0.6)
        )

        grid::grid.text(
          label = cell$label,
          x = grid::unit(cell$text_x, "npc"),
          y = grid::unit(cell$text_y, "npc"),
          just = c("left", "center"),
          gp = grid::gpar(col = self$textCol, cex = self$textCex, fontface = self$fontface)
        )
      }

      invisible()
    }
  )
)


# ---- SeqLegend
# This is the key change you requested: SeqLegend() builds on SeqBlankTrack(),
# wiping aesthetics first, then placing LegendKeys via SeqLegendElement.

SeqLegend <- function(...,
                      ncol = NULL, nrow = NULL,
                      colSpacing = 0, rowSpacing = 0,
                      hjust = 0, vjust = 0.5,
                      keyWidth = 0.04, keyGap = 0.006,
                      rowHeight = 0.12,
                      textCex = 0.9, textCol = "#1C1B1A",
                      fontface = 1,
                      height = NULL,
                      blankAesthetics = list()) {

  trk <- SeqBlankTrack(
    windows = .SeqUnitWindow(),     # <-- force ONE window for legend tracks
    height = height,
    aesthetics = blankAesthetics,
    singleWindow = TRUE
  )

  leg_el <- SeqLegendElement$new(
    ...,
    ncol = ncol, nrow = nrow,
    colSpacing = colSpacing, rowSpacing = rowSpacing,
    hjust = hjust, vjust = vjust,
    keyWidth = keyWidth, keyGap = keyGap,
    rowHeight = rowHeight,
    textCex = textCex, textCol = textCol,
    fontface = fontface
  )

  if (!is.null(trk$addElement) && is.function(trk$addElement)) {
    trk$addElement(leg_el)
  } else {
    trk$elements <- append(trk$elements, list(leg_el))
  }

  trk
}