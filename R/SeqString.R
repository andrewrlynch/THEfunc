# SeqString ----
#' SeqString: Strand-aware structural variant string links
#'
#' \code{SeqString} draws smooth cubic Bézier "string" curves linking
#' two genomic loci (breakpoints), typically representing SVs across
#' copy-number states.
#'
#' When \code{type="auto"}, curve type is inferred from strand:
#' \itemize{
#'   \item \code{+/+} → \code{"c"}
#'   \item \code{-/-} → \code{"c"}
#'   \item \code{+/-} → \code{"s"}
#'   \item \code{-/+} → \code{"s"}
#' }
#'
#' @param gr1 \code{GRanges}. First breakpoint set.
#' @param gr2 \code{GRanges}. Second breakpoint set.
#' @param t0,t1 Track indices (0 means "use arc track").
#' @param y0,y1 Numeric y anchor values in data units (e.g., copy number).
#' @param type \code{"auto"}, \code{"c"}, or \code{"s"} (vectorized/recycled).
#' @param curvature \code{"length"}, \code{"equal"}, or numeric (vectorized/recycled).
#' @param color Line color (vectorized/recycled).
#' @param width Line width (vectorized/recycled).
#' @param alpha Alpha (vectorized/recycled).
#' @param orientation Orientation (kept for compatibility; vectorized/recycled).
#'
#' @return An R6 \code{SeqString} object.
#' @export
SeqString <- R6::R6Class("SeqString",
                         inherit = SeqLink,
                         public = list(

                           type = NULL,
                           curvature = NULL,
                           aesthetics = list(),

                           #' @description
                           #' Create a new `SeqString` object.
                           #' @param gr1,gr2 GRange objects of start (`gr1`) and end (`gr2`) positions of SV junctions.
                           #' @param aesthetics A named list of aesthetics to override defaults (`shape`, `size`, `color`).
                           #' @return A new `SeqString` object.
                           initialize = function(gr1, gr2,
                                                 t0 = 0, t1 = 0,
                                                 y0 = 0, y1 = 0,
                                                 type = "auto",
                                                 curvature = "length",
                                                 color = "purple",
                                                 width = 1.5,
                                                 alpha = 0.8,
                                                 orientation = "*") {

                             super$initialize(
                               gr1 = gr1, gr2 = gr2, x1 = x1, x2 = x2,
                               t0 = t0, t1 = t1,
                               y0 = y0, y1 = y1,
                               color = color,
                               orientation = orientation
                             )

                             N <- length(self$gr1)

                             self$type      <- .recycled_to_N(type, N, "auto")
                             self$curvature <- .recycled_to_N(curvature, N, "length")
                             self$orientation <- .recycled_to_N(self$orientation, N, "*")

                             self$aesthetics <- list(
                               color = .recycled_to_N(color, N, "purple"),
                               width = .recycled_to_N(width, N, 1.5),
                               alpha = .recycled_to_N(alpha, N, 0.8)
                             )
                           },

                           #' @description
                           #' Prepare link coordinates in canvas space and infer strand-based curve types.
                           #' @param layout_all_tracks Layout metadata for all tracks/windows.
                           #' @param track_windows_list List of per-track \code{GRanges} windows.
                           #' @param arc_track_idx Track index used when \code{t0} or \code{t1} is 0.
                           prep = function(layout_all_tracks, track_windows_list, arc_track_idx) {

                             N <- length(self$gr1)
                             if (N == 0) return(invisible(NULL))

                             # resolve track indices (0 -> arc track)
                             t0 <- ifelse(self$t0 == 0, arc_track_idx, self$t0)
                             t1 <- ifelse(self$t1 == 0, arc_track_idx, self$t1)

                             ov1_all <- rep(NA_integer_, N)
                             ov2_all <- rep(NA_integer_, N)

                             # overlaps for left endpoints
                             for (tid in unique(t0)) {
                               idxs <- which(t0 == tid)
                               ov_matches <- GenomicRanges::findOverlaps(self$gr1[idxs], track_windows_list[[tid]], select = "all")
                               if (length(ov_matches) > 0) {
                                 matched_q <- idxs[S4Vectors::queryHits(ov_matches)]
                                 matched_s <- S4Vectors::subjectHits(ov_matches)
                                 ov1_all[matched_q] <- matched_s
                               }
                             }

                             # overlaps for right endpoints
                             for (tid in unique(t1)) {
                               idxs <- which(t1 == tid)
                               ov_matches <- GenomicRanges::findOverlaps(self$gr2[idxs], track_windows_list[[tid]], select = "all")
                               if (length(ov_matches) > 0) {
                                 matched_q <- idxs[S4Vectors::queryHits(ov_matches)]
                                 matched_s <- S4Vectors::subjectHits(ov_matches)
                                 ov2_all[matched_q] <- matched_s
                               }
                             }

                             valid <- !is.na(ov1_all) & !is.na(ov2_all)
                             if (!any(valid)) {
                               self$coordGrid <- NULL
                               return(invisible(NULL))
                             }

                             idx_valid <- which(valid)
                             nv <- length(idx_valid)

                             # normalize per-link params AFTER filtering
                             type_vec <- .recycled_to_N(self$type, N, "auto")[idx_valid]
                             type_vec <- tolower(as.character(type_vec))
                             auto_idx <- which(is.na(type_vec) | type_vec == "" | type_vec == "auto")
                             if (length(auto_idx) > 0) {
                               type_vec[auto_idx] <- stringTypeFromStrand(
                                 self$gr1[idx_valid][auto_idx],
                                 self$gr2[idx_valid][auto_idx],
                                 default = "c"
                               )
                             }
                             type_vec[!type_vec %in% c("c", "s")] <- "c"

                             curvature_vec <- .recycled_to_N(self$curvature, N, "length")[idx_valid]
                             orient_vec    <- .recycled_to_N(self$orientation, N, "*")[idx_valid]
                             col_vec       <- .recycled_to_N(self$aesthetics$color, N, "purple")[idx_valid]
                             lwd_vec       <- .recycled_to_N(self$aesthetics$width, N, 1.5)[idx_valid]
                             alpha_vec     <- .recycled_to_N(self$aesthetics$alpha, N, 0.8)[idx_valid]

                             x0s <- y0s <- x1s <- y1s <- numeric(nv)

                             for (k in seq_len(nv)) {
                               i <- idx_valid[k]

                               win1 <- ov1_all[i]
                               win2 <- ov2_all[i]

                               lay0 <- layout_all_tracks[[t0[i]]][[win1]]
                               lay1 <- layout_all_tracks[[t1[i]]][[win2]]

                               x0_gen <- IRanges::start(self$gr1[i])
                               x1_gen <- IRanges::start(self$gr2[i])

                               uv0 <- convertDataToGrid(x0_gen, self$y0[i], lay0$xscale, lay0$yscale)
                               uv1 <- convertDataToGrid(x1_gen, self$y1[i], lay1$xscale, lay1$yscale)

                               p0 <- gridToCanvas(uv0[1], uv0[2], lay0)
                               p1 <- gridToCanvas(uv1[1], uv1[2], lay1)

                               x0s[k] <- p0[1]; y0s[k] <- p0[2]
                               x1s[k] <- p1[1]; y1s[k] <- p1[2]
                             }

                             s1_vec <- as.character(GenomicRanges::strand(self$gr1[idx_valid]))
                             s2_vec <- as.character(GenomicRanges::strand(self$gr2[idx_valid]))


                             self$coordGrid <- data.frame(
                               x0 = x0s, y0 = y0s,
                               x1 = x1s, y1 = y1s,
                               strand1 = s1_vec,
                               strand2 = s2_vec,
                               type = type_vec,
                               curvature = curvature_vec,
                               orientation = orient_vec,
                               col = col_vec,
                               lwd = lwd_vec,
                               alpha = alpha_vec,
                               stringsAsFactors = FALSE
                             )

                             invisible(NULL)
                           },

                           #' @description
                           #' Draw all prepared string curves.
                           draw = function() {
                             df <- self$coordGrid
                             if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

                             for (i in seq_len(nrow(df))) {
                               drawSeqString(
                                 x0 = df$x0[i], y0 = df$y0[i],
                                 x1 = df$x1[i], y1 = df$y1[i],
                                 strand1 = df$strand1[i],
                                 strand2 = df$strand2[i],
                                 lwd = df$lwd[i],
                                 col = df$col[i],
                                 alpha = df$alpha[i]
                               )
                             }

                             invisible(NULL)
                           }
                         )
)