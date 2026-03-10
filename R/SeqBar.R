# SeqBar ----
#' SeqBar R6 Class
#'
#' @description
#' R6 class for drawing bar plots from genomic intervals in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' Each genomic interval is represented as a bar spanning its start–end
#' coordinates on the x-axis, with height determined by a y-column or a
#' default constant. Bars can be stacked by group if a grouping column is
#' provided. A default color palette is automatically assigned to groups
#' if no fill colors are specified.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100, 200), width = 50),
#'   score = c(2, 5, 3),
#'   group = c("A", "B", "A")
#' )
#' bars <- SeqBar(gr, yCol = "score", groupCol = "group")
#' bars$prep(layout_track = some_layout, track_windows = some_windows)
#' bars$draw()
#'
#' @export
SeqBar <- R6::R6Class("SeqBar",
                      inherit = SeqElement,
                      public = list(

                        #' @field gr A `GRanges` object containing genomic intervals.
                        gr = NULL,

                        #' @field yCol Optional column name in `gr` used for bar heights.
                        yCol = NULL,

                        #' @field groupCol Optional column name in `gr` defining groups for stacked bars.
                        groupCol = NULL,

                        #' @field groupLevels Optional character vector specifying factor levels for groups.
                        groupLevels = NULL,

                        #' @field y Numeric vector of bar heights (derived from `yCol` or constant).
                        y = NULL,

                        #' @field yStackedMax Maximum stacked y-value across groups, used for scaling.
                        yStackedMax = NULL,

                        #' @field group Factor defining group membership of each bar.
                        group = NULL,

                        #' @field aesthetics List of current aesthetics merged with defaults.
                        aesthetics = NULL,

                        #' @field coordCanvas List of data.frames containing transformed bar coordinates
                        #'   and fill colors for each window.
                        coordCanvas = NULL,

                        #' @field defaultAesthetics Default aesthetics for bar drawing:
                        #'   \code{fill = "grey60"}, \code{col = "#1C1B1A"}, \code{width = 0.8}, \code{lwd = 1}.
                        defaultAesthetics = list(
                          fill = "grey60",
                          col = "#1C1B1A",
                          width = 0.8,
                          lwd = 1
                        ),

                        #' @description
                        #' Create a new `SeqBar` object.
                        #'
                        #' @param gr A `GRanges` object with genomic intervals.
                        #' @param yCol Optional column name in `gr` for bar heights.
                        #' @param groupCol Optional column name in `gr` for grouping bars.
                        #' @param groupLevels Optional vector of group levels to enforce order.
                        #' @param aesthetics Optional list of aesthetic overrides.
                        #' @return A new `SeqBar` object.
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

                          self$y <- if (!is.null(yCol) && yCol %in% names(mcols(gr_use))) {
                            as.numeric(mcols(gr_use)[[yCol]])
                          } else {
                            rep(1, length(gr_use))
                          }

                          self$group <- if (!is.null(groupCol) && groupCol %in% names(mcols(gr_use))) {
                            as.character(mcols(gr_use)[[groupCol]])
                          } else {
                            rep("default", length(gr_use))
                          }

                          if (!is.null(groupLevels)) {
                            self$group <- factor(self$group, levels = groupLevels)
                          } else {
                            self$group <- factor(self$group)
                          }

                          # Auto-assign fill colors by group
                          if (is.null(self$aesthetics$fillPalette)) {
                            pal <- flexoki_palette(length(levels(self$group)))
                            names(pal) <- levels(self$group)
                            self$aesthetics$fillPalette <- pal
                          }

                          # Compute stacked maximum height
                          if (!is.null(self$groupCol)) {
                            xmid <- (start(self$gr) + end(self$gr)) / 2
                            group <- if (!is.null(groupLevels)) {
                              factor(as.character(mcols(self$gr)[[groupCol]]), levels = groupLevels)
                            } else {
                              as.factor(mcols(self$gr)[[groupCol]])
                            }

                            df <- data.frame(x = xmid, y = self$y, group = group)
                            df <- df[order(df$x, df$group), ]
                            y_cum <- tapply(df$y, df$x, sum)
                            self$yStackedMax <- max(y_cum, na.rm = TRUE)
                          } else {
                            self$yStackedMax <- max(self$y, na.rm = TRUE)
                          }

                          if (!is.null(aes)) {
                            resolved <- .resolve_aes(
                              data_mcols    = as.data.frame(S4Vectors::mcols(gr)),
                              aes_obj       = aes, scale_obj = scale,
                              n             = length(gr), default_color = "#1C1B1A"
                            )
                            # aes fill overrides fillPalette for grouped coloring
                            if (!is.null(resolved$fill)) {
                              self$aesthetics$fill <- resolved$fill
                            }
                            for (nm in setdiff(names(resolved), "fill")) self$aesthetics[[nm]] <- resolved[[nm]]
                          }
                        },

                        #' @description
                        #' Prepare bar coordinates in canvas space for each genomic window.
                        #'
                        #' @param layout_track A list of panel layout metadata for a track.
                        #' @param track_windows A `GRanges` object defining genomic windows.
                        prep = function(layout_track, track_windows) {
                          self$coordCanvas <- vector("list", length(track_windows))
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

                            start_pos <- start(gr_win)
                            end_pos   <- end(gr_win)
                            xmid      <- (start_pos + end_pos) / 2

                            df <- data.frame(
                              x = xmid,
                              x0 = start_pos,
                              x1 = end_pos,
                              y = y_win,
                              group = group_win,
                              stringsAsFactors = FALSE
                            )

                            clip <- clipToXscale(df$x0, df$x1, p$xscale)
                            if (length(clip$x0) == 0) {
                              self$coordCanvas[[w]] <- NULL
                              next
                            }

                            df <- df[clip$mask, ]
                            df$x0 <- clip$x0
                            df$x1 <- clip$x1

                            df$group <- factor(df$group, levels = levels(self$group))
                            df <- df[order(df$x, df$group), ]

                            if (!is.null(self$groupCol)) {
                              df$y0 <- 0
                              df$y1 <- 0
                              for (x in unique(df$x)) {
                                idx <- which(df$x == x)
                                heights <- df$y[idx]
                                df$y0[idx] <- cumsum(c(0, head(heights, -1)))
                                df$y1[idx] <- cumsum(heights)
                              }
                            } else {
                              df$y0 <- 0
                              df$y1 <- df$y
                            }

                            u0 <- (df$x0 - p$xscale[1]) / diff(p$xscale)
                            u1 <- (df$x1 - p$xscale[1]) / diff(p$xscale)
                            xleft  <- p$inner$x0 + u0 * (p$inner$x1 - p$inner$x0)
                            xright <- p$inner$x0 + u1 * (p$inner$x1 - p$inner$x0)

                            v0 <- (df$y0 - p$yscale[1]) / diff(p$yscale)
                            v1 <- (df$y1 - p$yscale[1]) / diff(p$yscale)
                            ybottom <- p$inner$y0 + v0 * (p$inner$y1 - p$inner$y0)
                            ytop <- p$inner$y0 + v1 * (p$inner$y1 - p$inner$y0)

                            if(!is.null(self$aesthetics$fill)){
                              fill_colors <- self$aesthetics$fill
                            } else {
                              fill_colors <- self$aesthetics$fillPalette[as.character(df$group)]
                            }

                            self$coordCanvas[[w]] <- data.frame(
                              x0 = xleft,
                              x1 = xright,
                              y0 = ybottom,
                              y1 = ytop,
                              fill = fill_colors,
                              stringsAsFactors = FALSE
                            )
                          }
                        },

                        #' @description
                        #' Draw bars onto the plotting canvas.
                        draw = function() {
                          if (is.null(self$coordCanvas)) return()
                          for (coords in self$coordCanvas) {
                            if (!is.data.frame(coords) || nrow(coords) == 0) next

                            grid.rect(
                              x = unit((coords$x0 + coords$x1) / 2, "npc"),
                              y = unit((coords$y0 + coords$y1) / 2, "npc"),
                              width = unit(abs(coords$x1 - coords$x0), "npc"),
                              height = unit(abs(coords$y1 - coords$y0), "npc"),
                              gp = gpar(
                                fill = coords$fill,
                                col = self$aesthetics$col,
                                lwd = self$aesthetics$lwd
                              )
                            )
                          }
                        }
                      )
)