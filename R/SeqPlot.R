# SeqPlot ----
#' SeqPlot R6 Class
#'
#' @description
#' The top-level container class in the SeqPlot framework. A `SeqPlot` manages
#' multiple genomic tracks, defines global windows, computes layout metadata,
#' and renders the complete multi-track plot including grid backgrounds, axes,
#' and sequence elements.
#'
#' @details
#' A `SeqPlot` object contains one or more `SeqTrack` objects, each of which
#' may contain `SeqElement` objects such as points, bars, lines, or links.
#' The `SeqPlot` class handles arranging these tracks into a grid, computing
#' coordinate transformations, applying aesthetics, and invoking draw routines.
#'
#' @examples
#' library(GenomicRanges)
#' win <- GRanges("chr1", IRanges(c(1, 1001), width = 500))
#' track1 <- SeqTrack(windows = win)
#' track2 <- SeqTrack(windows = win)
#' sp <- SeqPlot(tracks = list(track1, track2), windows = win)
#' sp$layoutGrid()
#' sp$drawGrid()
#' sp$drawAxes()
#' sp$drawElements()
#'
#' @export
SeqPlot <- R6Class("SeqPlot",
                   public = list(
                     #' @field tracks List of `SeqTrack` objects contained in the plot.
                     tracks = NULL,
                     #' @field windows Global genomic windows as a `GRanges` object.
                     windows = NULL,
                     #' @field layout Layout metadata produced by `$layoutGrid()`, containing
                     #' panel and track bounds.
                     layout = NULL,
                     #' @field aesthetics Named list of global aesthetics controlling plot-wide
                     #' appearance, merged with `defaultAesthetics`.
                     aesthetics = NULL,
                     #' @field defaultAesthetics Default aesthetics for track and window layout,
                     #' backgrounds, borders, axis lines, ticks, labels, and titles.
                     defaultAesthetics = list(
                       trackHeights = 1,
                       trackGaps = 0.01,
                       windowGaps = 0.01,
                       margins = list(top = 0.05, right = 0.05, bottom = 0.05, left = 0.05),
                       trackBackground = NA,
                       trackBorder = NA,
                       windowBackground = "whitesmoke",
                       windowBorder = "grey50",
                       xAxisLine = TRUE,
                       yAxisLine = TRUE,
                       xAxisBreakLines = FALSE,
                       yAxisBreakLines = FALSE,
                       xAxisTicks = TRUE,
                       yAxisTicks = TRUE,
                       xAxisLabels = TRUE,
                       yAxisLabels = TRUE,
                       xAxisLabelRotation = 0,
                       xAxisLabelVerticalJust = 1,
                       xAxisLabelHorizontalJust = 0.5,
                       xAxisTitle = TRUE,
                       yAxisTitle = TRUE,
                       yAxisTitleRotation = 0,
                       yAxisTitleVerticalJust = 0.5,
                       yAxisTitleHorizontalJust = 1,
                       yAxisPerWindow = FALSE,
                       axisColor = "#1C1B1A"
                     ),


                     #' @description
                     #' Create a new `SeqPlot` object.
                     #'
                     #' @param tracks List of `SeqTrack` objects.
                     #' @param windows Global `GRanges` windows. Defaults to `defaultGenomeWindows()`.
                     #' @param layout Optional layout object (normally produced by `$layoutGrid()`).
                     #' @param aesthetics Named list of aesthetics overriding `defaultAesthetics`.
                     #' @return A new `SeqPlot` object.
                     initialize = function(tracks = list(), windows = defaultGenomeWindows(), layout = NULL, aesthetics = list()) {
                       self$tracks <- tracks
                       self$windows <- windows
                       self$layout <- layout
                       self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)
                     },

                     #' @description
                     #' Compute the layout grid for all tracks and windows, assigning panel
                     #' coordinates, track heights, and axis scales.
                     #' @return Updates the `layout` field.
                     layoutGrid = function() {

                       .expand_limits <- function(r, expand) {
                         r <- range(r, na.rm = TRUE)
                         if (!all(is.finite(r))) r <- c(0, 1)
                         if (r[1] == r[2]) r <- r + c(-0.5, 0.5)  # avoid zero-span
                         span <- diff(r)
                         mul  <- if (length(expand) >= 1) expand[1] else 0
                         add  <- if (length(expand) >= 2) expand[2] else 0
                         pad  <- mul * span + add
                         list(data = r, expanded = c(r[1] - pad, r[2] + pad))
                       }

                       # --- Scale inheritance ---
                       for (i in seq_along(self$tracks)) {
                         track_i <- self$tracks[[i]]
                         # Inherit windows from SeqPlot if needed
                         if (is.null(track_i$windows) && is.null(track_i$scale_x)) {
                           if (is.null(self$windows)) {
                             stop("Global SeqPlot windows and at least one SeqTrack windows are NULL. All track windows must be set.")
                           }
                           track_i$windows  <- self$windows
                           track_i$scale_x  <- seq_scale_genomic(self$windows)
                         }
                         # If track has windows but no explicit scale_x, build one
                         if (!is.null(track_i$windows) && is.null(track_i$scale_x)) {
                           track_i$scale_x <- seq_scale_genomic(track_i$windows)
                         }
                       }

                       # --- Auto-infer scale_y from elements ---
                       for (i in seq_along(self$tracks)) {
                         track_i <- self$tracks[[i]]
                         if (is.null(track_i$scale_y)) {
                           for (elem in track_i$elements) {
                             inferred <- elem$.infer_scale_y()
                             if (!is.null(inferred)) {
                               track_i$scale_y <- inferred
                               break
                             }
                           }
                         }
                       }

                       nTracks <- length(self$tracks)
                       trackHeights <- if (length(self$aesthetics$trackHeights) == 1) rep(self$aesthetics$trackHeights, nTracks) else self$aesthetics$trackHeights

                       if (length(self$aesthetics$trackGaps) == 1) {
                         trackGaps <- rep(self$aesthetics$trackGaps, nTracks - 1)
                       } else if (length(self$aesthetics$trackGaps) == nTracks - 1) {
                         trackGaps <- self$aesthetics$trackGaps
                       } else {
                         stop("trackGaps must be a scalar or a vector of length (number of tracks - 1).")
                       }

                       if (length(self$aesthetics$windowGaps) == 1) {
                         windowGaps <- rep(self$aesthetics$windowGaps, nTracks)
                       } else if (length(self$aesthetics$windowGaps) == nTracks) {
                         windowGaps <- self$aesthetics$windowGaps
                       } else {
                         stop("windowGaps must be a scalar or a vector of length (number of tracks)")
                       }

                       margins <- self$aesthetics$margins
                       dataWidths <- lapply(self$tracks, function(t) width(t$windows))

                       # Calculate panel-level parameters
                       windowWidths <- lapply(self$tracks, function(t) {
                         win <- t$windows
                         raw_widths <- width(win)

                         # Check for per-window scale factor
                         if (!"scale" %in% names(mcols(win))) {
                           mcols(win)$scale <- rep(1e-6, length(win))  # default to Mb
                         }

                         effective_widths <- raw_widths * mcols(win)$scale
                         rel_widths <- effective_widths / sum(effective_widths)
                         return(rel_widths)
                       })

                       # xscales <- unlist(lapply(self$tracks, function(t) {
                       #   lapply(seq_along(t$windows), function(i) {
                       #     c(start(t$windows)[i], end(t$windows)[i])
                       #   })
                       # }), recursive = FALSE)

                       yscales <- list()

                       for (track_idx in seq_along(self$tracks)) {
                         track <- self$tracks[[track_idx]]
                         n_windows <- length(track$windows)
                         sy <- track$scale_y

                         # --- Check for genomic y_windows first (highest priority) ---
                         if (!is.null(track$y_windows) && track$uses_genomic_y) {
                           # Genomic y-axis from y_windows
                           n_y_win <- length(track$y_windows)
                           # If single y_window, apply to all x-windows; otherwise paired
                           if (n_y_win == n_windows) {
                             # Per-window y_windows
                             for (yw_idx in seq_len(n_y_win)) {
                               y_range <- c(start(track$y_windows)[yw_idx], end(track$y_windows)[yw_idx])
                               yscales <- c(yscales, list(y_range))
                             }
                           } else if (n_y_win == 1) {
                             # Single y_window applies to all x-windows
                             y_range <- c(start(track$y_windows)[1], end(track$y_windows)[1])
                             yscales <- c(yscales, rep(list(y_range), n_windows))
                           } else {
                             # Multiple y-windows apply to all x-windows (2D grid layout).
                             # y_sub_panels will handle per-window positioning; here just compute overall range.
                             y_overall <- c(min(GenomicRanges::start(track$y_windows)), 
                                           max(GenomicRanges::end(track$y_windows)))
                             yscales <- c(yscales, rep(list(y_overall), n_windows))
                           }
                           next
                         }

                         # --- Dispatch on explicit scale_y ---
                         if (!is.null(sy) && sy$type == "genomic") {
                           y_range <- c(min(GenomicRanges::start(sy$windows)),
                                        max(GenomicRanges::end(sy$windows)))
                           yscales <- c(yscales, rep(list(y_range), n_windows))
                           next

                         } else if (!is.null(sy) && sy$type == "discrete") {
                           n_lvl   <- length(sy$levels)
                           y_range <- c(0.5, n_lvl + 0.5)
                           yscales <- c(yscales, rep(list(y_range), n_windows))
                           next

                         } else if (!is.null(sy) && sy$type == "continuous" &&
                                    !is.null(sy$limits)) {
                           yscales <- c(yscales, rep(list(sy$limits), n_windows))
                           next
                         }

                         # --- Legacy path: auto-range from element data ---
                         disjoint <- isTRUE(track$aesthetics$disjointYScale)

                         if (!is.null(track$aesthetics$yAxisLimits)) {
                           # Single limit (apply to all windows)
                           if (is.numeric(track$aesthetics$yAxisLimits) &&
                               length(track$aesthetics$yAxisLimits) == 2) {
                             for (w in seq_len(n_windows)) {
                               yscales[[length(yscales) + 1]] <- track$aesthetics$yAxisLimits
                             }
                             next
                           }

                           # Vector of per-window limits
                           if (is.list(track$aesthetics$yAxisLimits) &&
                               length(track$aesthetics$yAxisLimits) == n_windows) {
                             for (w in seq_len(n_windows)) {
                               lim <- track$aesthetics$yAxisLimits[[w]]
                               yscales[[length(yscales) + 1]] <- if (length(lim) == 2) lim else c(0, 1)
                             }
                             next
                           }
                         } else if (disjoint) {
                           for (w in seq_len(n_windows)) {
                             win_gr <- track$windows[w]
                             y_vals <- c()

                             for (elem in track$elements) {
                               if (!is.null(elem$gr)) {
                                 ov <- findOverlaps(elem$gr, win_gr)
                                 if (length(ov) == 0) next

                                 idx <- S4Vectors::queryHits(ov)

                                 if (!is.null(elem$y))   y_vals <- c(y_vals, elem$y[idx])
                                 if (!is.null(elem$y0))  y_vals <- c(y_vals, elem$y0[idx])
                                 if (!is.null(elem$y1))  y_vals <- c(y_vals, elem$y1[idx])
                                 if (!is.null(elem$yStackedMax))  y_vals <- c(y_vals, elem$yStackedMax)
                               }

                               if (!is.null(elem$gr1) && !is.null(elem$gr2)) {
                                 ov1 <- findOverlaps(elem$gr1, win_gr)
                                 ov2 <- findOverlaps(elem$gr2, win_gr)
                                 idx <- unique(c(S4Vectors::queryHits(ov1), S4Vectors::queryHits(ov2)))

                                 if (!is.null(elem$y0))  y_vals <- c(y_vals, elem$y0[idx])
                                 if (!is.null(elem$y1))  y_vals <- c(y_vals, elem$y1[idx])
                                 if (!is.null(elem$height)) y_vals <- c(y_vals, elem$height[idx])
                               }
                             }

                             y_vals <- y_vals[is.finite(y_vals)]
                             if (length(y_vals) > 0) {
                               y_min <- min(y_vals)
                               y_max <- max(y_vals)
                               if (y_min == y_max) {
                                 y_min <- y_min - 1
                                 y_max <- y_max + 1
                               }
                               yscales[[length(yscales) + 1]] <- c(y_min, y_max)
                             } else {
                               yscales[[length(yscales) + 1]] <- c(0, 1)
                             }
                           }

                         } else {
                           win_gr <- track$windows
                           y_vals <- c()

                           for (elem in track$elements) {
                             if (!is.null(elem$gr)) {
                               ov <- GenomicRanges::findOverlaps(elem$gr, win_gr)
                               if (length(ov) == 0) next

                               idx <- S4Vectors::queryHits(ov)

                               if (!is.null(elem$y))   y_vals <- c(y_vals, elem$y[idx])
                               if (!is.null(elem$y0))  y_vals <- c(y_vals, elem$y0[idx])
                               if (!is.null(elem$y1))  y_vals <- c(y_vals, elem$y1[idx])
                               if (!is.null(elem$yStackedMax))  y_vals <- c(y_vals, elem$yStackedMax)
                             }

                             if (!is.null(elem$gr1) && !is.null(elem$gr2)) {
                               ov1 <- GenomicRanges::findOverlaps(elem$gr1, win_gr)
                               ov2 <- GenomicRanges::findOverlaps(elem$gr2, win_gr)
                               idx <- unique(c(S4Vectors::queryHits(ov1), S4Vectors::queryHits(ov2)))

                               if (!is.null(elem$y0))  y_vals <- c(y_vals, elem$y0[idx])
                               if (!is.null(elem$y1))  y_vals <- c(y_vals, elem$y1[idx])
                               if (!is.null(elem$height)) y_vals <- c(y_vals, elem$height[idx])
                             }
                           }

                           y_vals <- y_vals[is.finite(y_vals)]
                           if (length(y_vals) > 0) {
                             y_min <- min(y_vals)
                             y_max <- max(y_vals)
                             if (y_min == y_max) {
                               y_min <- y_min - 1
                               y_max <- y_max + 1
                             }
                           } else {
                             y_min <- 0
                             y_max <- 1
                           }

                           yscales <- c(yscales, rep(list(c(y_min, y_max)), n_windows))
                         }
                       }


                       availWidth <- 1 - (margins$left + margins$right)
                       availHeight <- 1 - (margins$top + margins$bottom)

                       total_gap_height <- sum(trackGaps)
                       availTrackHeight <- availHeight - total_gap_height
                       relHeights <- trackHeights / sum(trackHeights)
                       trackHeights_rel <- relHeights * availTrackHeight

                       y_top <- 1 - margins$top
                       grid_meta <- list()
                       panel_index <- 1

                       y_tops <- numeric(nTracks)
                       cursor <- 1 - margins$top

                       for (t in seq_len(nTracks)) {
                         if (t > 1) cursor <- cursor - trackGaps[t - 1]
                         y_tops[t] <- cursor
                         cursor <- cursor - trackHeights_rel[t]

                         h <- trackHeights_rel[t]
                         y_top_track <- y_tops[t]
                         y_bottom_track <- y_top_track - trackHeights_rel[t]


                         nWin <- length(windowWidths[[t]])
                         this_gap <- windowGaps[t]
                         gap_vec <- if (nWin == 1) numeric(0) else rep(this_gap, nWin - 1)
                         total_gap <- sum(gap_vec)

                         availWinWidth <- availWidth - total_gap
                         winFrac <- windowWidths[[t]]

                         expandX <- (self$tracks[[t]]$aesthetics$xExpand %||% self$aesthetics$xExpand) %||% c(0.05, 0)
                         expandY <- (self$tracks[[t]]$aesthetics$yExpand %||% self$aesthetics$yExpand) %||% c(0.05, 0)


                         grid_meta[[t]] <- vector("list", nWin)
                         x_left <- margins$left

                         # Pre-compute y_sub_panels for multiple y-windows
                         track <- self$tracks[[t]]
                         y_sub_panels_template <- NULL
                         if (!is.null(track$y_windows) && track$uses_genomic_y && length(track$y_windows) > 1) {
                           # Create y-sub-panels for the multiple y-windows
                           y_win_all <- track$y_windows
                           n_yw <- length(y_win_all)
                           if (!"scale" %in% names(mcols(y_win_all))) {
                             mcols(y_win_all)$scale <- rep(1e-6, n_yw)
                           }
                           y_win_gap <- track$aesthetics$yWindowGap %||% 0
                           total_y_gap <- y_win_gap * (n_yw - 1)
                           raw_h <- width(y_win_all) * mcols(y_win_all)$scale
                           avail_h <- (y_top_track - y_bottom_track) - total_y_gap
                           rel_h <- (raw_h / sum(raw_h)) * avail_h
                           y_sub_panels_template <- vector("list", n_yw)
                           y_cur <- y_bottom_track
                           for (yw in seq_len(n_yw)) {
                             sub_y0 <- y_cur
                             sub_y1 <- y_cur + rel_h[yw]
                             y_sub_panels_template[[yw]] <- list(
                               y0 = sub_y0, y1 = sub_y1,
                               yscale = c(start(y_win_all)[yw], end(y_win_all)[yw]),
                               yScaleFactor = mcols(y_win_all)$scale[[yw]],
                               y_seqname = as.character(seqnames(y_win_all)[yw])
                             )
                             y_cur <- sub_y1 + y_win_gap
                           }
                         }

                         for (w in seq_len(nWin)) {
                           w_width <- winFrac[w] * availWinWidth
                           x0_win <- x_left
                           x1_win <- x0_win + w_width

                           panel_full <- list(
                             x0 = x0_win,
                             x1 = x1_win,
                             y0 = y_bottom_track,
                             y1 = y_top_track
                           )

                           # --- NEW: compute per-window x data range and expand it
                           x_data <- c(start(track$windows)[w], end(track$windows)[w])
                           x_ex   <- .expand_limits(x_data, expandX)

                           # --- NEW: take the precomputed yscale for this panel_index and expand it
                           y_data <- yscales[[panel_index]]
                           y_ex   <- .expand_limits(y_data, expandY)

                           grid_meta[[t]][[w]] <- list(
                             track  = t,
                             window = w,
                             full   = panel_full,
                             inner  = panel_full,
                             # store both true ranges and expanded ranges
                             data_x = x_ex$data,
                             data_y = y_ex$data,
                             xscale = x_ex$expanded,
                             yscale = y_ex$expanded,
                             xScaleFactor = mcols(track$windows)$scale[[w]],
                             # scale type metadata
                             x_scale_type = if (!is.null(track$scale_x)) track$scale_x$type else "genomic",
                             y_scale_type = if (!is.null(track$y_windows) && track$uses_genomic_y) "genomic"
                                            else if (!is.null(track$scale_y)) track$scale_y$type
                                            else "continuous",
                             y_levels     = if (!is.null(track$scale_y) && track$scale_y$type == "discrete")
                                              track$scale_y$levels else NULL,
                             y_labels     = if (!is.null(track$scale_y) && track$scale_y$type == "discrete")
                                              (track$scale_y$labels %||% track$scale_y$levels) else NULL,
                             yScaleFactor = if (!is.null(track$y_windows) && track$uses_genomic_y) {
                                              # Get per-window scale factor from y_windows metadata
                                              mcols(track$y_windows)$scale[[min(w, length(track$y_windows))]]
                                            } else if (!is.null(track$scale_y) && track$scale_y$type == "genomic") {
                                              track$scale_y$scale_factor[min(w, length(track$scale_y$scale_factor))]
                                            } else NULL,
                             y_is_genomic = !is.null(track$y_windows) && track$uses_genomic_y,
                             y_sub_panels = y_sub_panels_template
                           )

                           if (w < nWin && length(gap_vec) >= w) {
                             x_left <- x_left + w_width + gap_vec[w]
                           }

                           panel_index <- panel_index + 1
                         }


                         y_top <- y_bottom_track - if (t == 1) 0 else trackGaps[t - 1]
                       }

                       trackBounds <- lapply(seq_along(grid_meta), function(t) {
                         panels <- grid_meta[[t]]
                         # extract the npc coords
                         x0s <- vapply(panels, function(p) p$full$x0, numeric(1))
                         x1s <- vapply(panels, function(p) p$full$x1, numeric(1))
                         y0s <- vapply(panels, function(p) p$full$y0, numeric(1))
                         y1s <- vapply(panels, function(p) p$full$y1, numeric(1))
                         list(
                           x0 = min(x0s), x1 = max(x1s),
                           y0 = min(y0s), y1 = max(y1s)
                         )
                       })

                       self$layout <- list(
                         panelBounds = grid_meta,
                         trackBounds = trackBounds
                       )

                       grid.newpage()
                       pushViewport(viewport(name = "root"))
                       popViewport()
                     },

                     #' @description
                     #' Draw the grid backgrounds, track areas, window panels, and borders.
                     #' @return Renders the grid to the graphics device.
                     drawGrid = function() {
                       stopifnot(is.list(self$layout),
                                 !is.null(self$layout$panelBounds),
                                 !is.null(self$layout$trackBounds))

                       panelBounds      <- self$layout$panelBounds
                       trackBounds <- self$layout$trackBounds
                       nTracks     <- length(panelBounds)

                       for (t in seq_len(nTracks)) {

                         trackAesthetics = modifyList(self$aesthetics, self$tracks[[t]]$aesthetics)

                         # Draw the full‐track background/border using trackBounds[t]:
                         tb <- trackBounds[[t]]
                         grid.rect(
                           x      = unit(tb$x0, "npc"),
                           y      = unit(tb$y0, "npc"),
                           width  = unit(tb$x1 - tb$x0, "npc"),
                           height = unit(tb$y1 - tb$y0, "npc"),
                           just   = c("left", "bottom"),
                           gp     = gpar(
                             fill = trackAesthetics$trackBackground,
                             col = trackAesthetics$trackBorder,
                             lwd = 0.5
                           )
                         )

                         for (win in panelBounds[[t]]) {
                           p <- win$full
                           xscale <- win$xscale
                           yscale <- win$yscale

                           if (!is.null(win$y_sub_panels)) {
                             # Multiple y-windows: draw sub-panel backgrounds and borders
                             for (sub in win$y_sub_panels) {
                               grid.rect(
                                 x = unit(p$x0, "npc"), y = unit(sub$y0, "npc"),
                                 width = unit(p$x1 - p$x0, "npc"), height = unit(sub$y1 - sub$y0, "npc"),
                                 just = c("left", "bottom"),
                                 gp = gpar(fill = trackAesthetics$windowBackground, col = NA, lwd = 0.5)
                               )
                             }
                             for (sub in win$y_sub_panels) {
                               grid.rect(
                                 x = unit(p$x0, "npc"), y = unit(sub$y0, "npc"),
                                 width = unit(p$x1 - p$x0, "npc"), height = unit(sub$y1 - sub$y0, "npc"),
                                 just = c("left", "bottom"),
                                 gp = gpar(fill = NA, col = trackAesthetics$windowBorder, lwd = 1)
                               )
                             }
                           } else {
                             # Draw window panel background inside the track
                             grid.rect(
                               x      = unit(p$x0, "npc"),
                               y      = unit(p$y0, "npc"),
                               width  = unit(p$x1 - p$x0, "npc"),
                               height = unit(p$y1 - p$y0, "npc"),
                               just   = c("left", "bottom"),
                               gp     = gpar(
                                 fill = trackAesthetics$windowBackground,
                                 col = NA,
                                 lwd = 0.5
                               )
                             )

                             # Draw window panel borders inside the track
                             grid.rect(
                               x      = unit(p$x0, "npc"),
                               y      = unit(p$y0, "npc"),
                               width  = unit(p$x1 - p$x0, "npc"),
                               height = unit(p$y1 - p$y0, "npc"),
                               just   = c("left", "bottom"),
                               gp     = gpar(
                                 fill = NA,
                                 col = trackAesthetics$windowBorder,
                                 lwd = 1
                               )
                             )
                           }

                           # Draw window x breaks over the background
                           if (isTRUE(trackAesthetics$xAxisBreakLines)) {
                             xbreaks <- xscale  # ← replace with custom logic later
                             xgrid   <- (xbreaks - xscale[1]) / diff(xscale)
                             for (x in xgrid) {
                               grid.lines(
                                 x = unit(c(win$full$x0 + x * (win$full$x1 - win$full$x0)), "npc"),
                                 y = unit(c(win$full$y0, win$full$y1), "npc"),
                                 gp = gpar(col = "grey90", lwd = 0.5)
                               )
                             }
                           }

                           # Draw window y breaks over the background
                           if (isTRUE(trackAesthetics$yAxisBreakLines)) {
                             ybreaks <- yscale
                             ygrid   <- (ybreaks - yscale[1]) / diff(yscale)
                             for (y in ygrid) {
                               grid.lines(
                                 x = unit(c(win$full$x0, win$full$x1), "npc"),
                                 y = unit(c(win$full$y0 + y * (win$full$y1 - win$full$y0)), "npc"),
                                 gp = gpar(col = "grey90", lwd = 0.5)
                               )
                             }
                           }

                           # Draw window panel borders inside the track
                           # Only draw outer border when no y_sub_panels exist;
                           # individual sub-panel borders drawn above in the y_sub_panels loop
                           if (is.null(win$y_sub_panels)) {
                             grid.rect(
                               x      = unit(p$x0, "npc"),
                               y      = unit(p$y0, "npc"),
                               width  = unit(p$x1 - p$x0, "npc"),
                               height = unit(p$y1 - p$y0, "npc"),
                               just   = c("left", "bottom"),
                               gp     = gpar(
                                 fill = NA,
                                 col = trackAesthetics$windowBorder,
                                 lwd = 1
                               )
                             )
                           }
                         }
                       }

                       #popViewport()
                     },

                     #' @description
                     #' Draw x- and y-axes for all tracks and windows, including ticks, labels,
                     #' and axis titles.
                     #' @return Renders axes to the graphics device.
                     drawAxes = function() {

                       # Helper functions
                       # Null-coalescing convenience
                       `%||%` <- function(a, b) if (!is.null(a)) a else b

                       # Format unit label from scale factor
                       .axis_unit_label <- function(sf) {
                         if (abs(sf - 1e-6) < 1e-9) {
                           "Mb"
                         } else if (abs(sf - 1e-3) < 1e-9) {
                           "kb"
                         } else if (abs(sf - 1) < 1e-9) {
                           "bp"
                         } else {
                           paste0("×", signif(sf, 2))
                         }
                       }

                       # Format numbers with big marks and adaptive rounding
                       .axis_format_num <- function(x, digits = NULL) {
                         if (is.null(digits)) {
                           rng <- diff(range(na.omit(x)))
                           digits <- if (is.finite(rng)) {
                             if (rng >= 100) 0 else if (rng >= 10) 1 else 2
                           } else 0
                         }
                         format(round(x, digits), big.mark = ",", scientific = FALSE, trim = TRUE)
                       }

                       # Core axis metadata generator
                       make_axis_meta <- function(plot_range,
                                                  data_range = NULL,     # true genomic/data range
                                                  manual_breaks = NULL,  # numeric or NULL
                                                  n = 5,                 # target # of pretty breaks
                                                  cap = c("full", "capped", "exact", "ticks"),
                                                  labels = NULL,
                                                  minor_breaks = NULL    # NULL, integer, or numeric
                       ) {
                         cap <- match.arg(cap)
                         pr <- range(plot_range, na.rm = TRUE)  # already expanded
                         dr <- if (!is.null(data_range)) range(data_range, na.rm = TRUE) else pr

                         # --- major breaks ---

                         br <- if (!is.null(manual_breaks)) {
                           as.numeric(manual_breaks)
                         } else {
                           # Use the true data range (if supplied) for stable breaks across windows
                           # and fall back to the expanded plot range otherwise.
                           br_range <- if (!is.null(data_range)) dr else pr
                           br_fun   <- scales::pretty_breaks(n = n)
                           br_fun(br_range)
                         }

                         # Always include 0 when it lies within the plotted range
                         if (0 >= pr[1] && 0 <= pr[2]) {
                           br <- c(0, br)
                         }

                         # Filter to the plotted range (with tiny tolerance) and de-duplicate
                         tol <- diff(pr) * 1e-9
                         br <- br[is.finite(br) & br >= (pr[1] - tol) & br <= (pr[2] + tol)]
                         br <- sort(unique(br))


                         # --- axis line range ---
                         ar <- switch(cap,
                                      full   = pr,
                                      capped = if (length(br)) range(br) else pr,
                                      exact = dr,
                                      ticks  = NULL)

                         # --- major labels ---
                         lab <- if (!is.null(labels)) labels else as.character(br)

                         # --- minor breaks ---
                         mbr <- NULL
                         if (!is.null(minor_breaks)) {
                           if (is.numeric(minor_breaks) && length(minor_breaks) > 1) {
                             # explicit numeric vector
                             mbr <- minor_breaks
                           } else if (is.numeric(minor_breaks) && length(minor_breaks) == 1) {
                             # evenly spaced subdivisions
                             nb <- as.integer(minor_breaks)
                             if (nb > 1 && length(br) > 1) {
                               mids <- unlist(lapply(seq_along(br)[-length(br)], function(i) {
                                 seq(br[i], br[i+1], length.out = nb + 1)[-c(1, nb+1)]
                               }))
                               mbr <- mids
                             }
                           }
                           mbr <- mbr[is.finite(mbr) & mbr >= pr[1] & mbr <= pr[2]]
                         }

                         list(
                           data_range   = dr,   # true data range (unexpanded)
                           expand_range = pr,   # plotting range (expanded, from layoutGrid)
                           breaks       = br,
                           labels       = lab,
                           minor_breaks = mbr,
                           axis_range   = ar
                         )
                       }



                       panelBounds      <- self$layout$panelBounds
                       trackBounds <- self$layout$trackBounds
                       nTracks     <- length(panelBounds)

                       defaults <- list(
                         # expansion
                         xExpand = c(0.025, 0),
                         yExpand = c(0.05, 0),

                         # major breaks
                         xAxisBreaks = NULL, yAxisBreaks = NULL,
                         xAxisNBreaks = 5,   yAxisNBreaks = 4,

                         # minor breaks
                         xMinorBreaks = NULL,   # NULL (none), integer (# subdivisions), or numeric vector
                         yMinorBreaks = NULL,

                         # capping
                         xAxisCap = "capped",
                         yAxisCap = "capped",

                         # styling
                         axisColor        = "#1C1B1A",
                         tickLengthNPC    = 0.005,
                         xTitleOffsetNPC  = 0.06,
                         yTitleOffsetNPC  = 0.05,
                         xLabelOffsetNPC  = 0.015,
                         yLabelOffsetNPC  = 0.010,
                         xAxisLabelSize   = 0.6,
                         yAxisLabelSize   = 0.6,

                         xAxisLabelRotation = 0,
                         xAxisLabelHorizontalJust = 0.5,
                         xAxisLabelVerticalJust   = 1
                       )

                       for (t in seq_len(nTracks)) {
                         trackAesthetics = modifyList(
                           modifyList(defaults, self$aesthetics %||% list()),
                           self$tracks[[t]]$aesthetics %||% list()
                           )

                         # Draw window axes over the panel borders
                         for (win in panelBounds[[t]]) {
                           xscale <- win$xscale
                           yscale <- win$yscale
                           p <- win$full

                           # axis metadata
                           x_meta <- make_axis_meta(
                             plot_range   = win$xscale,
                             data_range   = win$data_x,
                             manual_breaks= trackAesthetics$xAxisBreaks,
                             n            = trackAesthetics$xAxisNBreaks,
                             cap          = trackAesthetics$xAxisCap,
                             minor_breaks = trackAesthetics$xMinorBreaks
                           )

                           y_meta <- make_axis_meta(
                             plot_range   = win$yscale,
                             data_range   = win$data_y,
                             manual_breaks= trackAesthetics$yAxisBreaks,
                             n            = trackAesthetics$yAxisNBreaks,
                             cap          = trackAesthetics$yAxisCap,
                             minor_breaks = trackAesthetics$yMinorBreaks
                           )

                           # --- Y-axis formatting dispatch ---
                           y_scale_type <- win$y_scale_type %||% "continuous"

                           if (y_scale_type == "genomic") {
                             sf_y       <- win$yScaleFactor %||% 1e-6
                             unit_lbl_y <- .axis_unit_label(sf_y)
                             scaled_y   <- y_meta$breaks * sf_y
                             ylabels    <- .axis_format_num(scaled_y)
                             if (length(ylabels) > 0)
                               ylabels[length(ylabels)] <- paste0(ylabels[length(ylabels)], " ", unit_lbl_y)
                             y_meta$labels <- ylabels

                           } else if (y_scale_type == "discrete") {
                             disc_n         <- length(win$y_levels)
                             y_meta$breaks  <- seq_len(disc_n)
                             y_meta$labels  <- win$y_labels %||% win$y_levels
                             y_meta$minor_breaks <- NULL
                             y_meta$axis_range   <- c(0.5, disc_n + 0.5)

                           } else {
                             # continuous — legacy yAxisGenomicLabels flag still honored
                             if (isTRUE(trackAesthetics$yAxisGenomicLabels)) {
                               sf_y       <- win$yScaleFactor %||% 1e-6
                               unit_lbl_y <- .axis_unit_label(sf_y)
                               scaled_y   <- y_meta$breaks * sf_y
                               ylabels    <- .axis_format_num(scaled_y)
                               if (length(ylabels) > 0)
                                 ylabels[length(ylabels)] <- paste0(ylabels[length(ylabels)], " ", unit_lbl_y)
                               y_meta$labels <- ylabels
                             }
                           }

                           map_x_npc <- function(val) {
                             (val - x_meta$expand_range[1]) / diff(x_meta$expand_range) * (p$x1 - p$x0) + p$x0
                           }
                           map_y_npc <- function(val) {
                             (val - y_meta$expand_range[1]) / diff(y_meta$expand_range) * (p$y1 - p$y0) + p$y0
                           }

                           # x‐axis along bottom of this track
                           if (isTRUE(trackAesthetics$xAxisLine) && !is.null(x_meta$axis_range)) {
                             x0_line <- map_x_npc(x_meta$axis_range[1])
                             x1_line <- map_x_npc(x_meta$axis_range[2])
                             grid.lines(
                               x = unit(c(x0_line, x1_line), "npc"),
                               y = unit(c(p$y0, p$y0), "npc"),
                               gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                             )
                           }

                           # ticks + labels on X
                           if (isTRUE(trackAesthetics$xAxisTicks) && length(x_meta$breaks)) {
                             x_scale_type <- win$x_scale_type %||% "genomic"

                             if (x_scale_type == "genomic") {
                               sf <- if (!is.null(win$xScaleFactor)) win$xScaleFactor else 1e-6
                               unit_label <- .axis_unit_label(sf)
                               scaled_vals <- x_meta$breaks * sf
                               xlabels <- .axis_format_num(scaled_vals)
                               if (length(xlabels) > 0) {
                                 xlabels[length(xlabels)] <- paste0(xlabels[length(xlabels)], " ", unit_label)
                               }
                             } else {
                               xlabels <- .axis_format_num(x_meta$breaks)
                             }

                             # minor ticks (drawn once, outside major tick loop)
                             if (!is.null(x_meta$minor_breaks) && length(x_meta$minor_breaks)) {
                               for (mb in x_meta$minor_breaks) {
                                 xpos_mb <- map_x_npc(mb)
                                 grid.lines(
                                   x = unit(c(xpos_mb, xpos_mb), "npc"),
                                   y = unit(c(p$y0, p$y0 - trackAesthetics$tickLengthNPC * 0.6), "npc"),
                                   gp = gpar(col = trackAesthetics$axisColor, lwd = 0.5)
                                 )
                               }
                             }

                             for (i in seq_along(x_meta$breaks)) {
                               xpos <- map_x_npc(x_meta$breaks[i])

                               # major tick
                               grid.lines(
                                 x = unit(c(xpos, xpos), "npc"),
                                 y = unit(c(p$y0, p$y0 - trackAesthetics$tickLengthNPC), "npc"),
                                 gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                               )

                               # label
                               if (isTRUE(trackAesthetics$xAxisLabels)) {
                                 grid.text(
                                   label = xlabels[i],
                                   x = unit(xpos, "npc"),
                                   y = unit(p$y0 - trackAesthetics$xLabelOffsetNPC, "npc"),
                                   just  = "top",
                                   rot   = trackAesthetics$xAxisLabelRotation,
                                   hjust = trackAesthetics$xAxisLabelHorizontalJust,
                                   vjust = trackAesthetics$xAxisLabelVerticalJust,
                                   gp = gpar(cex = trackAesthetics$xAxisLabelSize)
                                 )
                               }
                             }
                           }

                           # X axis title
                           if (isTRUE(trackAesthetics$xAxisTitle)) {
                             x_label <- if (!is.null(trackAesthetics$xAxisTitleText) &&
                                            !isTRUE(trackAesthetics$xAxisTitleText)) {
                               trackAesthetics$xAxisTitleText
                             } else if ((win$x_scale_type %||% "genomic") == "genomic") {
                               gsub("^chr", "", as.character(seqnames(self$tracks[[t]]$windows[win$window])))
                             } else {
                               ""
                             }
                             if (nzchar(x_label)) {
                               grid.text(
                                 label = x_label,
                                 x = unit((p$x0 + p$x1) / 2, "npc"),
                                 y = unit(p$y0 - trackAesthetics$xTitleOffsetNPC, "npc"),
                                 just = "top",
                                 gp = gpar(cex = 0.6, fontface = "bold")
                               )
                             }
                           }

                           # --- Y AXIS (left of panel or first window only) ---
                           draw_y_here <- isTRUE(trackAesthetics$yAxisLine) &&
                             (isTRUE(trackAesthetics$yAxisPerWindow) || win$window == 1)

                           if (!is.null(win$y_sub_panels) && draw_y_here) {
                             # Multiple y-windows: draw axis per sub-panel
                             for (sub in win$y_sub_panels) {
                               sf_y <- sub$yScaleFactor %||% 1e-6
                               unit_label_y <- .axis_unit_label(sf_y)
                               sub_y_meta <- make_axis_meta(
                                 plot_range   = sub$yscale,
                                 data_range   = sub$yscale,
                                 manual_breaks = trackAesthetics$yAxisBreaks,
                                 n            = trackAesthetics$yAxisNBreaks,
                                 cap          = trackAesthetics$yAxisCap,
                                 minor_breaks = trackAesthetics$yMinorBreaks
                               )
                               map_sub_y <- function(val) {
                                 (val - sub$yscale[1]) / diff(sub$yscale) * (sub$y1 - sub$y0) + sub$y0
                               }
                               if (isTRUE(trackAesthetics$yAxisLine) && !is.null(sub_y_meta$axis_range)) {
                                 grid.lines(
                                   x = unit(c(p$x0, p$x0), "npc"),
                                   y = unit(c(map_sub_y(sub_y_meta$axis_range[1]),
                                              map_sub_y(sub_y_meta$axis_range[2])), "npc"),
                                   gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                                 )
                               }
                               if (isTRUE(trackAesthetics$yAxisTicks) && length(sub_y_meta$breaks)) {
                                 scaled_y_vals <- sub_y_meta$breaks * sf_y
                                 ylabels <- .axis_format_num(scaled_y_vals)
                                 if (length(ylabels) > 0) {
                                   ylabels[length(ylabels)] <- paste0(ylabels[length(ylabels)], " ", unit_label_y)
                                 }
                                 for (i in seq_along(sub_y_meta$breaks)) {
                                   ypos <- map_sub_y(sub_y_meta$breaks[i])
                                   grid.lines(
                                     x = unit(c(p$x0, p$x0 - trackAesthetics$tickLengthNPC), "npc"),
                                     y = unit(c(ypos, ypos), "npc"),
                                     gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                                   )
                                   if (isTRUE(trackAesthetics$yAxisLabels)) {
                                     grid.text(
                                       label = ylabels[i],
                                       x = unit(p$x0 - trackAesthetics$yLabelOffsetNPC, "npc"),
                                       y = unit(ypos, "npc"),
                                       just = "right",
                                       gp = gpar(cex = trackAesthetics$yAxisLabelSize)
                                     )
                                   }
                                 }
                               }
                             }
                           } else {
                             if (draw_y_here && !is.null(y_meta$axis_range)) {
                               y0_line <- map_y_npc(y_meta$axis_range[1])
                               y1_line <- map_y_npc(y_meta$axis_range[2])
                               grid.lines(
                                 x = unit(c(p$x0, p$x0), "npc"),
                                 y = unit(c(y0_line, y1_line), "npc"),
                                 gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                               )
                             }

                             if (isTRUE(trackAesthetics$yAxisTicks) &&
                                 (isTRUE(trackAesthetics$yAxisPerWindow) || win$window == 1) &&
                                 length(y_meta$breaks)) {

                               # minor ticks (drawn once, outside major tick loop)
                               if (!is.null(y_meta$minor_breaks) && length(y_meta$minor_breaks)) {
                                 for (mb in y_meta$minor_breaks) {
                                   ypos_mb <- map_y_npc(mb)
                                   grid.lines(
                                     x = unit(c(p$x0, p$x0 - trackAesthetics$tickLengthNPC * 0.6), "npc"),
                                     y = unit(c(ypos_mb, ypos_mb), "npc"),
                                     gp = gpar(col = trackAesthetics$axisColor, lwd = 0.5)
                                   )
                                 }
                               }

                               for (i in seq_along(y_meta$breaks)) {
                                 ypos <- map_y_npc(y_meta$breaks[i])

                                 # major tick
                                 grid.lines(
                                   x = unit(c(p$x0, p$x0 - trackAesthetics$tickLengthNPC), "npc"),
                                   y = unit(c(ypos, ypos), "npc"),
                                   gp = gpar(col = trackAesthetics$axisColor, lwd = 1)
                                 )

                                 # label
                                 if (isTRUE(trackAesthetics$yAxisLabels)) {
                                   y_lab_text <- if (!is.null(y_meta$labels) && i <= length(y_meta$labels))
                                     y_meta$labels[i]
                                   else
                                     .axis_format_num(y_meta$breaks[i])
                                   grid.text(
                                     label = y_lab_text,
                                     x = unit(p$x0 - trackAesthetics$yLabelOffsetNPC, "npc"),
                                     y = unit(ypos, "npc"),
                                     just = "right",
                                     gp = gpar(cex = trackAesthetics$yAxisLabelSize)
                                   )
                                 }
                               }
                             }
                           }

                           # Y axis title
                           if (isTRUE(trackAesthetics$yAxisTitle) && win$window == 1) {
                             y_label <- if (!is.null(trackAesthetics$yAxisTitleText)) {
                               trackAesthetics$yAxisTitleText
                             } else if (length(self$tracks[[t]]$elements) > 0 &&
                                        !is.null(self$tracks[[t]]$elements[[1]]$yCol)) {
                               self$tracks[[t]]$elements[[1]]$yCol
                             } else {
                               ""
                             }

                             grid.text(
                               label = y_label,
                               x = unit(p$x0 - trackAesthetics$yTitleOffsetNPC, "npc"),
                               y = unit((p$y0 + p$y1) / 2, "npc"),
                               just = "center",
                               rot  = trackAesthetics$yAxisTitleRotation %||% 90,
                               hjust= trackAesthetics$yAxisTitleHorizontalJust %||% 0.5,
                               vjust= trackAesthetics$yAxisTitleVerticalJust   %||% 0.5,
                               gp = gpar(cex = trackAesthetics$yAxisLabelSize, fontface = "bold")
                             )
                           }
                         }
                       }
                     },

                     #' @description
                     #' Draw all elements from every track, including features (`SeqElement`)
                     #' and links (`SeqLink`).
                     #' @return Renders elements to the graphics device.
                     drawElements = function() {
                       # Plotting links first
                       track_windows_list <- lapply(self$tracks, function(trk) trk$windows)

                       for (track_idx in seq_along(self$tracks)) {
                         track <- self$tracks[[track_idx]]
                         for (elem in track$elements) {
                           if (is(elem, "SeqLink")) {
                             elem$prep(layout_all_tracks = self$layout$panelBounds,
                                       track_windows_list = track_windows_list,
                                       arc_track_idx = track_idx)
                             elem$draw()
                           }
                         }
                       }

                       # Plotting features
                       for (track_idx in seq_along(self$tracks)) {
                         track <- self$tracks[[track_idx]]
                         layout_track <- self$layout$panelBounds[[track_idx]]
                         track_windows <- track$windows

                         for (elem in track$elements) {
                           if (!is(elem, "SeqLink")){
                             if (is.element("prep", names(elem))) elem$prep(layout_track, track_windows)
                             if (is.element("draw", names(elem))) elem$draw()
                           }
                         }
                       }
                     },

                     #' @description
                     #' Append a `SeqTrack` to this plot.
                     #' @param track A `SeqTrack` object.
                     #' @return The `SeqPlot` object (invisibly).
                     addTrack = function(track) {
                       self$tracks <- append(self$tracks, list(track))
                       invisible(self)
                     },

                     #' @description
                     #' Prepares and draws all components of the SeqPlot.
                     #' @return Renders elements to the graphics device.
                     plot = function() {
                       self$layoutGrid()
                       self$drawGrid()
                       self$drawAxes()
                       self$drawElements()
                     }

                   ))
