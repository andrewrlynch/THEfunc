# SeqArea ----
#' SeqArea R6 Class
#'
#' @description
#' R6 class for drawing filled area plots from genomic intervals in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' Each genomic interval is represented as an area under a line, with y-values
#' taken from a metadata column or set to a constant. Areas can be grouped and
#' stacked, with each group assigned its own fill color. Missing group–x
#' combinations are padded with zero values to ensure continuous stacked shapes.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100, 200), width = 50),
#'   score = c(2, 5, 3),
#'   group = c("A", "B", "A")
#' )
#' area <- SeqArea(gr, yCol = "score", groupCol = "group")
#' area$prep(layout_track = some_layout, track_windows = some_windows)
#' area$draw()
#'
#' @export
SeqArea <- R6::R6Class("SeqArea",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic intervals.
                         gr = NULL,

                         #' @field yCol Optional column name in `gr` used for y-values.
                         yCol = NULL,

                         #' @field groupCol Optional column name in `gr` defining groups for stacked areas.
                         groupCol = NULL,

                         #' @field groupLevels Optional character vector specifying factor levels for groups.
                         groupLevels = NULL,

                         #' @field y Numeric vector of y-values (from `yCol` or constant).
                         y = NULL,

                         #' @field group Factor defining group membership of each interval.
                         group = NULL,

                         #' @field aesthetics List of current aesthetics merged with defaults.
                         aesthetics = NULL,

                         #' @field coordCanvas List of polygon coordinate lists for each group in
                         #'   each window, storing transformed x, y, and fill values.
                         coordCanvas = NULL,

                         #' @field yStackedMax Maximum stacked y-value across groups, used for scaling.
                         yStackedMax = NULL,

                         #' @field defaultAesthetics Default aesthetics for area drawing:
                         #'   \code{fill = "grey60"}, \code{col = "black"}, \code{alpha = 1}, \code{lwd = 0.5}.
                         defaultAesthetics = list(
                           fill = "grey60",
                           col = "black",
                           alpha = 1,
                           lwd = 0.5
                         ),

                         #' @description
                         #' Create a new `SeqArea` object.
                         #'
                         #' @param gr A `GRanges` object containing genomic intervals.
                         #' @param yCol Optional column name in `gr` for y-values.
                         #' @param groupCol Optional column name in `gr` for grouping areas.
                         #' @param groupLevels Optional vector of group levels to enforce order.
                         #' @param aesthetics Optional list of aesthetic overrides.
                         #' @return A new `SeqArea` object.
                         initialize = function(gr = NULL, x = NULL, yCol = NULL, groupCol = NULL,
                                               groupLevels = NULL, aesthetics = list(),
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
                           self$groupCol <- groupCol
                           self$groupLevels <- groupLevels
                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)

                           self$y <- if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
                             as.numeric(mcols(gr)[[yCol]])
                           } else {
                             rep(1, length(gr))
                           }

                           self$group <- if (!is.null(groupCol) && groupCol %in% names(mcols(gr))) {
                             as.character(mcols(gr)[[groupCol]])
                           } else {
                             rep("default", length(gr))
                           }

                           if (!is.null(groupLevels)) {
                             self$group <- factor(self$group, levels = groupLevels)
                           } else {
                             self$group <- factor(self$group)
                           }

                           # Compute stacked maximum height for scaling
                           xmid <- (start(gr) + end(gr)) / 2
                           df <- data.frame(x = xmid, y = self$y, group = self$group)
                           y_totals <- tapply(df$y, df$x, sum)
                           self$yStackedMax <- max(y_totals, na.rm = TRUE)

                           # Auto-assign fill colors by group
                           if (is.null(self$aesthetics$fillPalette)) {
                             pal <- flexoki_palette(length(levels(self$group)))
                             names(pal) <- levels(self$group)
                             self$aesthetics$fillPalette <- pal
                           }

                           if (!is.null(aes)) {
                             resolved <- .resolve_aes(
                               data_mcols    = as.data.frame(S4Vectors::mcols(gr)),
                               aes_obj       = aes, scale_obj = scale,
                               n             = length(gr), default_color = "#1C1B1A"
                             )
                             if (!is.null(resolved$fill)) self$aesthetics$fill <- resolved$fill
                             for (nm in setdiff(names(resolved), "fill")) self$aesthetics[[nm]] <- resolved[[nm]]
                           }
                         },

                         #' @description
                         #' Prepare stacked area polygons in canvas space for each genomic window.
                         #'
                         #' @param layout_track A list of panel layout metadata for a track.
                         #' @param track_windows A `GRanges` object defining genomic windows.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- list()

                           ov <- findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- queryHits(ov)
                           sh <- subjectHits(ov)
                           gr_sub <- self$gr[qh]
                           y_sub <- self$y[qh]
                           group_sub <- self$group[qh]

                           for (w in unique(sh)) {
                             p <- layout_track[[w]]
                             mask <- sh == w
                             idx <- qh[mask]
                             if (length(idx) == 0) next

                             gr_win <- self$gr[idx]
                             y_win <- y_sub[mask]
                             group_win <- group_sub[mask]
                             x <- (start(gr_win) + end(gr_win)) / 2

                             df <- data.frame(
                               x = x,
                               y = y_win,
                               group = group_win,
                               stringsAsFactors = FALSE
                             )
                             df$group <- factor(df$group, levels = levels(self$group))

                             # Pad missing group–x combinations
                             df <- tidyr::complete(df, x, group = levels(self$group), fill = list(y = 0))
                             df <- df[order(df$x, df$group), ]

                             # Compute stacked y0/y1 per group at each x
                             df$y0 <- NA_real_
                             df$y1 <- NA_real_
                             for (xval in unique(df$x)) {
                               rows <- which(df$x == xval)
                               running_y <- 0
                               for (i in rows) {
                                 df$y0[i] <- running_y
                                 running_y <- running_y + df$y[i]
                                 df$y1[i] <- running_y
                               }
                             }

                             # Convert to canvas coordinates
                             u <- (df$x - p$xscale[1]) / diff(p$xscale)
                             x_abs <- p$inner$x0 + u * (p$inner$x1 - p$inner$x0)

                             v0 <- (df$y0 - p$yscale[1]) / diff(p$yscale)
                             v1 <- (df$y1 - p$yscale[1]) / diff(p$yscale)
                             y0_abs <- p$inner$y0 + v0 * (p$inner$y1 - p$inner$y0)
                             y1_abs <- p$inner$y0 + v1 * (p$inner$y1 - p$inner$y0)

                             fill_colors <- self$aesthetics$fillPalette[as.character(df$group)]

                             groups <- split(
                               data.frame(x = x_abs, y0 = y0_abs, y1 = y1_abs, fill = fill_colors),
                               df$group
                             )

                             for (g in groups) {
                               g <- g[order(g$x), ]
                               x_poly <- c(g$x, rev(g$x))
                               y_poly <- c(g$y0, rev(g$y1))
                               fill_poly <- c(g$fill, rev(g$fill))

                               self$coordCanvas[[length(self$coordCanvas) + 1]] <- list(
                                 x = x_poly,
                                 y = y_poly,
                                 fill = fill_poly
                               )
                             }
                           }
                         },

                         #' @description
                         #' Draw stacked areas onto the plotting canvas.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()
                           for (poly in self$coordCanvas) {
                             grid.polygon(
                               x = unit(poly$x, "npc"),
                               y = unit(poly$y, "npc"),
                               gp = gpar(
                                 fill = poly$fill,
                                 col = self$aesthetics$col,
                                 lwd = self$aesthetics$lwd,
                                 alpha = self$aesthetics$alpha
                               )
                             )
                           }
                         }
                       )
)