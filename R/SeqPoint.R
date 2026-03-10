# SeqPoint ----
#' SeqPoint R6 Class
#'
#' @description
#' An R6 class for plotting genomic points, such as SNPs or single-base
#' features, on a SeqPlot track. Each genomic range is drawn as a point,
#' with optional y-axis values from metadata.
#'
#' @param gr GRange object of point coordinates.
#' @param y A numeric vector of y-values for each point. Defaults to 0.5 if no `yCol` is provided.
#' @param yCol Optional character string naming a metadata column in `gr` that will be used for y-values.
#' @param color Optional character string naming a metadata column in `gr` that will be used for color values.
#'
#' @export
SeqPoint <- R6::R6Class("SeqPoint",
                        inherit = SeqElement,
                        public = list(

                          #' @field gr A `GRanges` object containing the genomic positions
                          #'   of the points and any associated metadata.
                          gr = NULL,

                          #' @field y A numeric vector of y-values for each point. Defaults to
                          #'   0.5 if no `yCol` is provided.
                          y = NULL,

                          #' @field yCol Optional character string naming a metadata column in
                          #'   `gr` that will be used for y-values.
                          yCol = NULL,

                          #' @field color Optional character string naming a metadata column in `gr` that will be used for color values.
                          color = NULL,

                          #' @field coordOriginal A `GRanges` object storing the original,
                          #'   unmodified genomic coordinates.
                          coordOriginal = NULL,

                          #' @field coordCanvas A list of per-panel coordinate matrices
                          #'   (x, y in canvas units) produced by `prep()`.
                          coordCanvas = NULL,

                          #' @field coordIndex A list of per-panel coordinate indices.
                          coordIndex = NULL,

                          #' @field aesthetics A list of plotting aesthetics for the points,
                          #'   merged from user input and defaults.
                          aesthetics = NULL,

                          #' @field defaultAesthetics Default aesthetics for points:
                          #'   - `shape`: plotting symbol (default 16)
                          #'   - `size`: point size scaling factor
                          #'   - `color`: point color
                          defaultAesthetics = list(
                            shape = 16,
                            size = 0.1
                            #color = "#1C1B1A"
                          ),

                          #' @description
                          #' Create a new `SeqPoint` object.
                          #' @param gr A `GRanges` object of point features.
                          #' @param yCol Optional column in `gr` metadata to use as y-values.
                          #' @param aesthetics A named list of aesthetics to override defaults
                          #'   (`shape`, `size`, `color`).
                          #' @return A new `SeqPoint` object.
                          #' @examples
                          #' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:5, width = 1))
                          #' pt <- SeqPoint(gr)
                          #' pt$prep(layout_track, track_windows)
                          #' pt$draw()
                          initialize = function(gr = NULL, x = NULL, yCol = NULL, color = NULL,
                                                aesthetics = list(),
                                                aes = NULL, scale = NULL) {
                            # Accept either x (preferred) or legacy gr
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

                            if (!is.null(yCol) && yCol %in% names(mcols(gr_use))) {
                              self$y <- as.numeric(mcols(gr_use)[[yCol]])
                            } else {
                              self$y <- rep(0.5, length(gr_use))
                            }

                            n <- length(gr_use)

                            if (!is.null(color) && is.character(color) && length(color) == 1 &&
                                color %in% names(mcols(gr_use))) {
                              # Mapping: column of gr
                              color_vals <- as.character(mcols(gr_use)[[color]])
                            } else if (!is.null(color) && length(color) == n) {
                              # Mapping: explicit vector of colors
                              color_vals <- as.character(color)
                            } else {
                              # No mapping provided → fill with NA for now
                              color_vals <- rep(NA_character_, n)
                            }

                            # Merge with defaults and user aesthetics
                            aes_list <- modifyList(self$defaultAesthetics, aesthetics)

                            # Replace NA colors with fixed color if provided
                            if ("color" %in% names(aes_list)) {
                              idx_na <- is.na(color_vals)
                              color_vals[idx_na] <- aes_list$color
                            }

                            # Save resolved aesthetics
                            self$aesthetics <- aes_list
                            self$aesthetics$color <- color_vals  # always a vector now

                            # aes() mapping overrides (higher priority than color=)
                            if (!is.null(aes)) {
                              resolved <- .resolve_aes(
                                data_mcols    = as.data.frame(S4Vectors::mcols(gr)),
                                aes_obj       = aes,
                                scale_obj     = scale,
                                n             = n,
                                default_color = "#1C1B1A"
                              )
                              for (nm in names(resolved)) self$aesthetics[[nm]] <- resolved[[nm]]
                            }
                          },

                          #' @description
                          #' Prepare the point coordinates for plotting in canvas space,
                          #' transforming genomic ranges to per-panel grid coordinates.
                          #' @param layout_track A list of panel metadata for the current track
                          #'   from `SeqPlot$layoutGrid()`.
                          #' @param track_windows A `GRanges` object defining the genomic
                          #'   windows for the track.
                          prep = function(layout_track, track_windows) {
                            self$coordCanvas <- vector("list", length(track_windows))

                            ov <- GenomicRanges::findOverlaps(self$gr, track_windows)
                            if (length(ov) == 0) return(invisible())

                            qh <- S4Vectors::queryHits(ov)
                            sh <- S4Vectors::subjectHits(ov)

                            x <- start(self$gr)[qh]
                            y <- self$y[qh]

                            for (w in unique(sh)) {
                              mask <- sh == w
                              if (sum(mask) == 0) next

                              x_sub <- x[mask]
                              y_sub <- y[mask]
                              panel_meta <- layout_track[[w]]

                              u <- (x_sub - panel_meta$xscale[1]) / diff(panel_meta$xscale)
                              v <- (y_sub - panel_meta$yscale[1]) / diff(panel_meta$yscale)
                              u <- pmax(pmin(u, 1), 0)
                              v <- pmax(pmin(v, 1), 0)

                              x_canvas <- panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0)
                              y_canvas <- panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)

                              self$coordCanvas[[w]] <- cbind(x_canvas, y_canvas)
                              self$coordIndex[[w]]  <- qh[sh == w]
                            }

                            invisible()
                          },

                          #' @description
                          #' Draw the points using grid graphics, applying aesthetics for
                          #' shape, color, and size.
                          draw = function() {
                            if (is.null(self$coordCanvas)) return()
                            for (w in seq_along(self$coordCanvas)) {
                              coords <- self$coordCanvas[[w]]
                              idx <- self$coordIndex[[w]]
                              col_vals <- self$aesthetics$color[idx]

                              if (is.null(coords)) next
                              grid.points(
                                x = unit(coords[, 1], "npc"),
                                y = unit(coords[, 2], "npc"),
                                pch = self$aesthetics$shape,
                                gp = gpar(col = self$aesthetics$color[idx], cex = self$aesthetics$size)
                              )
                            }
                          }
                        )
)