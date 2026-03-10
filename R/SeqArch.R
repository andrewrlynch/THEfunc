# SeqArch ----
#' SeqArch R6 Class
#'
#' @description
#' R6 class for drawing arch-style links between genomic loci in the SeqPlot
#' framework. Inherits from [SeqLink].
#'
#' @details
#' `SeqArch` connects pairs of genomic intervals (`gr1`, `gr2`) with curved
#' arches that can span within or across tracks. Heights can be fixed or taken
#' from a metadata column (`yCol`). Stubs are supported for links where only
#' one endpoint lies within a visible window. Aesthetics control the colors
#' and widths of stems, arches, and stubs.
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 50), width = 1))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), width = 1))
#' arch <- SeqArch(gr1, gr2, t0 = 1, t1 = 1, height = 0.5)
#' # arch$prep(layout_all_tracks = some_layout, track_windows_list = some_windows, arc_track_idx = 1)
#' # arch$draw()
#'
#' @export
SeqArch <- R6::R6Class("SeqArch",
                       inherit = SeqLink,
                       public = list(

                         #' @field yCol Optional column name in `gr1` to set per-link arch height.
                         yCol = NULL,

                         #' @field height Numeric vector of arch heights (constant or from `yCol`).
                         height = NULL,

                         #' @field curve Character vector controlling curvature type, e.g. `"length"`
                         #'   (scale with span) or `"equal"` (uniform).
                         curve = NULL,

                         #' @field left_only Integer indices of links with only the left end visible.
                         left_only = NULL,

                         #' @field right_only Integer indices of links with only the right end visible.
                         right_only = NULL,

                         #' @field stubs List of stub objects representing partial arches.
                         stubs = NULL,

                         #' @field aesthetics List of aesthetics controlling stem/arch appearance:
                         #'   \itemize{
                         #'     \item stemColor — color of link stems
                         #'     \item arcColor — color of arches
                         #'     \item stemWidth — line width of stems
                         #'     \item arcWidth — line width of arches
                         #'     \item stubAngle — angle of stub hooks (degrees)
                         #'     \item stubLength — length of stub hooks (NPC units)
                         #'     \item plotStubs — logical; whether to draw stubs
                         #'   }
                         aesthetics = list(),

                         #' @description
                         #' Create a new `SeqArch` object.
                         #'
                         #' @param gr1 Deprecated; use `x1` instead.
                         #' @param gr2 Deprecated; use `x2` instead.
                         #' @param x1 Preferred: `GRanges` object giving start positions of arches.
                         #' @param x2 Preferred: `GRanges` object giving end positions of arches.
                         #' @param t0 Integer or vector giving track indices for start loci (default: `0`).
                         #' @param t1 Integer or vector giving track indices for end loci (default: `0`).
                         #' @param y0 Numeric or vector giving vertical start positions (default: `0`).
                         #' @param y1 Numeric or vector giving vertical end positions (default: `0`).
                         #' @param yCol Optional column name in `gr1` used to determine arch height.
                         #' @param orientation Character or vector specifying link orientation (default: `"*"`).
                         #' @param height Numeric scalar default height if `yCol` not used.
                         #' @param curve Character or vector specifying curvature mode (`"length"`, `"equal"`, or numeric).
                         #' @param aesthetics List of aesthetics for stems, arches, and stubs.
                         #' @return A new `SeqArch` object.
                         initialize = function(gr1 = NULL, gr2 = NULL, x1 = NULL, x2 = NULL, t0 = 0, t1 = 0,
                                               y0 = 0, y1 = 0, yCol = NULL, orientation = "*",
                                               height = 1, curve = "length", aesthetics = list()) {

                           super$initialize(gr1 = gr1, gr2 = gr2, x1 = x1, x2 = x2, t0 = t0, t1 = t1, y0 = y0, y1 = y1,
                                            color = aesthetics$color %||% "black", orientation = orientation)

                           if (!is.null(yCol) && yCol %in% names(mcols(gr1))) {
                             self$height <- as.numeric(mcols(gr1)[[yCol]])
                           } else {
                             self$height <- rep(height, length(gr1))
                           }

                           self$curve <- rep(curve, length(gr1))

                           self$aesthetics <- list(
                             stemColor = rep(aesthetics$stemColor %||% "black", length(gr1)),
                             arcColor  = rep(aesthetics$arcColor  %||% "black", length(gr1)),
                             stemWidth = rep(aesthetics$stemWidth %||% 1, length(gr1)),
                             arcWidth  = rep(aesthetics$arcWidth  %||% 1, length(gr1)),
                             stubAngle  = aesthetics$stubAngle  %||% 45,
                             stubLength = aesthetics$stubLength %||% 0.02,
                             plotStubs  = aesthetics$plotStubs  %||% TRUE
                           )
                         },

                         #' @description
                         #' Prepare arch coordinates and stubs for drawing.
                         #'
                         #' @param layout_all_tracks A nested list of panel layouts for all tracks.
                         #' @param track_windows_list A list of `GRanges` objects defining genomic windows per track.
                         #' @param arc_track_idx Integer index of the track where arches should be drawn.
                         prep = function(layout_all_tracks, track_windows_list, arc_track_idx) {
                           N <- length(self$gr1)
                           if (N == 0) return(NULL)

                           t0 <- ifelse(self$t0 == 0, arc_track_idx, self$t0)
                           t1 <- ifelse(self$t1 == 0, arc_track_idx, self$t1)

                           ov1_all <- rep(NA_integer_, N)
                           ov2_all <- rep(NA_integer_, N)

                           for (tid in unique(t0)) {
                             idxs <- which(t0 == tid)
                             ov_matches <- findOverlaps(self$gr1[idxs], track_windows_list[[tid]], select = "all")
                             if (length(ov_matches) > 0) {
                               matched_q <- idxs[queryHits(ov_matches)]
                               matched_s <- subjectHits(ov_matches)
                               ov1_all[matched_q] <- matched_s
                             }
                           }

                           for (tid in unique(t1)) {
                             idxs <- which(t1 == tid)
                             ov_matches <- findOverlaps(self$gr2[idxs], track_windows_list[[tid]], select = "all")
                             if (length(ov_matches) > 0) {
                               matched_q <- idxs[queryHits(ov_matches)]
                               matched_s <- subjectHits(ov_matches)
                               ov2_all[matched_q] <- matched_s
                             }
                           }

                           left_only  <- which(!is.na(ov1_all) &  is.na(ov2_all))
                           right_only <- which(is.na(ov1_all)  & !is.na(ov2_all))
                           self$left_only  <- left_only
                           self$right_only <- right_only

                           t0i  <- ifelse(self$t0==0, arc_track_idx, self$t0)
                           t1i  <- ifelse(self$t1==0, arc_track_idx, self$t1)
                           ov1  <- ov1_all
                           ov2  <- ov2_all

                           if (!isTRUE(self$aesthetics$plotStubs)) {
                             self$stubs <- NULL
                           } else {
                             self$stubs <- list()
                             for (i in c(left_only, right_only)) {
                               use1 <- i %in% left_only
                               tid_local <- if (use1) t0i[i] else t1i[i]
                               win_local <- if (use1) ov1[i]  else ov2[i]
                               pm_local  <- layout_all_tracks[[tid_local]][[win_local]]

                               gr         <- if (use1) self$gr1[i] else self$gr2[i]
                               yval_local <- if (use1) self$y0[i]  else self$y1[i]
                               uv_base    <- convertDataToGrid(start(gr), yval_local, pm_local$xscale, pm_local$yscale)
                               p_base     <- gridToCanvas(uv_base[1], uv_base[2], pm_local)

                               pm_arc     <- layout_all_tracks[[arc_track_idx]][[win_local]]
                               uv_top     <- convertDataToGrid(start(gr), self$height[i], pm_arc$xscale, pm_arc$yscale)
                               p_top      <- gridToCanvas(uv_top[1], uv_top[2], pm_arc)

                               dir_value <- if (i %in% left_only) +1 else -1
                               partner_chr <- if (use1) as.character(seqnames(self$gr2[i])) else as.character(seqnames(self$gr1[i]))

                               self$stubs[[length(self$stubs) + 1]] <- list(
                                 x     = p_base[1],
                                 y0    = p_base[2],
                                 y1    = p_top[2],
                                 dir   = dir_value,
                                 partner = partner_chr,
                                 color = self$aesthetics$stemColor[i],
                                 width = self$aesthetics$stemWidth[i]
                               )
                             }
                           }

                           has_any <- !is.na(ov1_all) | !is.na(ov2_all)
                           if (!any(has_any)) return()

                           valid <- !is.na(ov1_all) & !is.na(ov2_all)
                           if (!any(valid)) return()

                           df <- vector("list", sum(valid))
                           j <- 1
                           for (i in which(valid)) {
                             win1 <- ov1_all[i]
                             win2 <- ov2_all[i]

                             layout_x0 <- layout_all_tracks[[t0[i]]][[win1]]
                             layout_x1 <- layout_all_tracks[[t1[i]]][[win2]]
                             layout_y0 <- layout_all_tracks[[t0[i]]][[win1]]
                             layout_y1 <- layout_all_tracks[[t1[i]]][[win2]]

                             x0_gen <- start(self$gr1[i])
                             x1_gen <- start(self$gr2[i])

                             uv0 <- convertDataToGrid(x0_gen, self$y0[i], layout_x0$xscale, layout_y0$yscale)
                             uv1 <- convertDataToGrid(x1_gen, self$y1[i], layout_x1$xscale, layout_y1$yscale)

                             p0 <- gridToCanvas(uv0[1], uv0[2], layout_y0)
                             p1 <- gridToCanvas(uv1[1], uv1[2], layout_y1)

                             layout_arc1 <- layout_all_tracks[[arc_track_idx]][[win1]]
                             layout_arc2 <- layout_all_tracks[[arc_track_idx]][[win2]]

                             uvh0 <- convertDataToGrid(x0_gen, self$height[i], layout_arc1$xscale, layout_arc1$yscale)
                             uvh1 <- convertDataToGrid(x1_gen, self$height[i], layout_arc2$xscale, layout_arc2$yscale)

                             top0 <- gridToCanvas(uvh0[1], uvh0[2], layout_arc1)[2]
                             top1 <- gridToCanvas(uvh1[1], uvh1[2], layout_arc2)[2]

                             df[[j]] <- data.frame(
                               x0          = p0[1],   y0       = p0[2],
                               x1          = p1[1],   y1       = p1[2],
                               top0        = top0,    top1     = top1,
                               orientation = self$orientation[i],
                               curve       = self$curve[i],
                               stemColor   = self$aesthetics$stemColor[i],
                               arcColor    = self$aesthetics$arcColor[i],
                               stemWidth   = self$aesthetics$stemWidth[i],
                               arcWidth    = self$aesthetics$arcWidth[i],
                               stringsAsFactors = FALSE
                             )
                             j <- j + 1
                           }
                           self$coordGrid <- do.call(rbind, df)
                         },

                         #' @description
                         #' Draw arches and stubs onto the plotting canvas.
                         draw = function() {
                           df <- self$coordGrid
                           if (!is.null(df)) {
                             for (i in seq_len(nrow(df))) {
                               drawSeqArch(
                                 x0 = df$x0[i], y0 = df$y0[i],
                                 x1 = df$x1[i], y1 = df$y1[i],
                                 top0 = df$top0[i], top1 = df$top1[i],
                                 orientation = df$orientation[i],
                                 curve = df$curve[i],
                                 stemWidth = df$stemWidth[i],
                                 stemColor = df$stemColor[i],
                                 arcWidth = df$arcWidth[i],
                                 arcColor = df$arcColor[i]
                               )
                             }
                           }

                           if (isTRUE(self$aesthetics$plotStubs) && !is.null(self$stubs) && length(self$stubs) > 0) {
                             angle_deg <- self$aesthetics$stubAngle %||% 45
                             theta     <- angle_deg * pi / 180
                             L_npc     <- self$aesthetics$stubLength %||% 0.02
                             dx        <- L_npc * cos(theta)
                             dy        <- L_npc * sin(theta)

                             xs   <- sapply(self$stubs, `[[`, "x")
                             y0s  <- sapply(self$stubs, `[[`, "y0")
                             y1s  <- sapply(self$stubs, `[[`, "y1")
                             dirs <- sapply(self$stubs, `[[`, "dir")
                             cols <- sapply(self$stubs, `[[`, "color")
                             wds  <- sapply(self$stubs, `[[`, "width")

                             grid.segments(
                               x0 = unit(xs,  "npc"),
                               y0 = unit(y0s, "npc"),
                               x1 = unit(xs,  "npc"),
                               y1 = unit(y1s, "npc"),
                               gp = gpar(
                                 col = scales::alpha(sapply(self$stubs, `[[`, "color"), 0.5),
                                 lwd = sapply(self$stubs, `[[`, "width")
                               )
                             )

                             grid.segments(
                               x0    = unit(xs,            "npc"),
                               y0    = unit(y1s,           "npc"),
                               x1    = unit(xs + dirs * dx, "npc"),
                               y1    = unit(y1s +            dy, "npc"),
                               gp    = gpar(col = scales::alpha(cols, 0.5), lwd = wds),
                               arrow = arrow(type = "open", length = unit(1, "mm"))
                             )

                             partners <- gsub("chr", "", sapply(self$stubs, `[[`, "partner"))
                             ys     <- sapply(self$stubs, `[[`, "y1")
                             xs     <- sapply(self$stubs, `[[`, "x")

                             grid.text(
                               label = partners,
                               x     = unit(xs, "npc"),
                               y     = unit(ys + 0.01, "npc"),
                               just  = c("center", "bottom"),
                               gp    = gpar(col = cols, fontsize = 7)
                             )
                           }
                         }
                       )
)