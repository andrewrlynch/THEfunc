
# SeqRecon ----
#' SeqRecon R6 Class
#'
#' @description
#' R6 class for plotting **ReCon-style** structural variant arches
#' (inversions, duplications/deletions, and translocations).
#' Inherits from [SeqArch].
#'
#' @details
#' `SeqRecon` extends `SeqArch` by automatically classifying structural
#' variants into three tiers:
#'
#' * **Inversion** (`+/+`, `-/-`) → head-to-head (HH) or tail-to-tail (TT).
#' * **Duplication/Deletion** (`-/+`, `+/-`) → tandem duplication or deletion.
#' * **Translocation** (different chromosomes).
#'
#' Each class is mapped to a fixed vertical tier and colored consistently.
#' Tier guide lines and labels are drawn automatically.
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 50), width = 1), strand = c("+", "-"))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), width = 1), strand = c("+", "-"))
#' recon <- SeqRecon(gr1, gr2)
#' # recon$prep(layout_all_tracks = some_layout, track_windows_list = some_windows, arc_track_idx = 1)
#' # recon$draw()
#'
#' @export
SeqRecon <- R6::R6Class("SeqRecon",
                        inherit = SeqArch,
                        public = list(

                          #' @field last_arc_track Cached layout of the last arc track used in `prep()`.
                          last_arc_track = NULL,

                          #' @field col_h2h Color for head-to-head inversions (`+/+`).
                          col_h2h = NULL,

                          #' @field col_t2t Color for tail-to-tail inversions (`-/-`).
                          col_t2t = NULL,

                          #' @field col_dup Color for tandem duplications (`-/+`).
                          col_dup = NULL,

                          #' @field col_del Color for deletions (`+/-`).
                          col_del = NULL,

                          #' @field col_trans Color for translocations (across chromosomes).
                          col_trans = NULL,

                          #' @field drawClasses Character vector of class names to draw
                          #'   (default: `c("Inversion", "Dup/Del", "Translocation")`).
                          drawClasses = c("Inversion", "Dup/Del", "Translocation"),

                          #' @field arc_track_idx Integer index of the track where arcs are drawn.
                          arc_track_idx = NULL,

                          #' @field layout_all_tracks Cached layout for all tracks used in `prep()`.
                          layout_all_tracks = NULL,

                          #' @field tierMultipliers Numeric vector assigning each variant to a fixed tier.
                          tierMultipliers = NULL,

                          #' @description
                          #' Create a new `SeqRecon` object.
                          #'
                          #' @param gr1 Deprecated start loci; use `x1` instead.
                          #' @param gr2 Deprecated end loci; use `x2` instead.
                          #' @param x1 Preferred: `GRanges` object with start loci.
                          #' @param x2 Preferred: `GRanges` object with end loci.
                          #' @param yCol Optional column name for custom arch heights.
                          #' @param y0,y1 Numeric start and end y-values (default: `0`).
                          #' @param t0,t1 Track indices for start and end loci (default: `0` → current track).
                          #' @param orientation Character vector of orientations (default: `"*"`).
                          #' @param curve Character curvature mode (`"length"`, `"equal"`, or numeric).
                          #' @param drawClasses Which class tiers to render.
                          #' @param aesthetics List of aesthetic overrides (colors, widths, etc.).
                          #' @return A new `SeqRecon` object.
                          initialize = function(gr1 = NULL, gr2 = NULL, x1 = NULL, x2 = NULL,
                                                yCol = NULL,
                                                y0 = 0, y1 = 0,
                                                t0 = 0, t1 = 0,
                                                orientation = "*",
                                                curve = "length",
                                                drawClasses = c("Inversion", "Dup/Del", "Translocation"),
                                                aesthetics = list()) {

                            self$col_h2h   <- aesthetics$h2hColor   %||% flexoki_palette(9)[3]
                            self$col_t2t   <- aesthetics$t2tColor   %||% flexoki_palette(9)[4]
                            self$col_dup   <- aesthetics$dupColor   %||% flexoki_palette(9)[1]
                            self$col_del   <- aesthetics$delColor   %||% flexoki_palette(9)[2]
                            self$col_trans <- aesthetics$transColor %||% flexoki_palette(9)[9]

                            super$initialize(
                              gr1        = gr1,
                              gr2        = gr2,
                              x1         = x1,
                              x2         = x2,
                              t0         = t0,
                              t1         = t1,
                              y0         = 0,
                              y1         = 0,
                              yCol       = yCol,
                              orientation= orientation,
                              height     = 1,
                              curve      = curve,
                              aesthetics = aesthetics
                            )

                            self$aesthetics <- modifyList(self$aesthetics, list(drawClasses = self$drawClasses))
                          },

                          #' @description
                          #' Prepare classification and layout for ReCon-style drawing.
                          #'
                          #' @param layout_all_tracks List of layout panels for all tracks.
                          #' @param track_windows_list List of genomic windows per track.
                          #' @param arc_track_idx Index of the arc track in which arches are drawn.
                          prep = function(layout_all_tracks, track_windows_list, arc_track_idx) {
                            self$last_arc_track <- layout_all_tracks[[arc_track_idx]]
                            self$coordGrid <- NULL

                            N    <- length(self$gr1)
                            seq1 <- as.character(seqnames(self$gr1))
                            seq2 <- as.character(seqnames(self$gr2))

                            s1   <- ifelse(as.character(strand(self$gr1)) %in% c("+","-"),
                                           as.character(strand(self$gr1)), "+")
                            s2   <- ifelse(as.character(strand(self$gr2)) %in% c("+","-"),
                                           as.character(strand(self$gr2)), "+")
                            code <- paste0(s1, "/", s2)

                            tier  <- numeric(N)
                            cols  <- character(N)
                            ori   <- character(N)

                            t0 <- ifelse(self$t0 == 0, arc_track_idx, self$t0)
                            t1 <- ifelse(self$t1 == 0, arc_track_idx, self$t1)

                            for (i in seq_len(N)) {
                              if (seq1[i] != seq2[i]) {
                                tier[i] <- 1
                                cols[i] <- self$col_trans
                                ori[i]  <- code[i]
                              } else if (code[i] %in% c("+/+", "-/-")) {
                                tier[i] <- 0
                                cols[i] <- if (code[i] == "+/+") self$col_h2h else self$col_t2t
                                ori[i]  <- ifelse(code[i] == "+/+", "+", "-")
                              } else if (code[i] %in% c("-/+", "+/-")) {
                                tier[i] <- 0.5
                                cols[i] <- if (code[i] == "-/+") self$col_dup else self$col_del
                                ori[i]  <- ifelse(code[i] == "-/+", "+", "-")
                              } else {
                                tier[i] <- 1
                                cols[i] <- self$col_trans
                                ori[i]  <- "+"
                              }
                            }

                            self$tierMultipliers <- tier
                            self$orientation <- ori
                            self$height      <- tier
                            self$aesthetics$arcColor  <- cols
                            self$aesthetics$stemColor <- cols
                            self$layout_all_tracks <- layout_all_tracks
                            self$arc_track_idx     <- arc_track_idx

                            super$prep(layout_all_tracks, track_windows_list, arc_track_idx)
                          },

                          #' @description
                          #' Draw ReCon-style arches, tiers, and labels.
                          draw = function() {
                            if (is.null(self$last_arc_track)) return()

                            panels <- self$last_arc_track
                            x0s <- vapply(panels, function(pm) pm$inner$x0, numeric(1))
                            x1s <- vapply(panels, function(pm) pm$inner$x1, numeric(1))
                            y0s <- vapply(panels, function(pm) pm$inner$y0, numeric(1))
                            y1s <- vapply(panels, function(pm) pm$inner$y1, numeric(1))
                            tb_x0 <- min(x0s); tb_x1 <- max(x1s)
                            tb_y0 <- min(y0s); tb_y1 <- max(y1s)
                            ysc   <- panels[[1]]$yscale

                            drawClasses = self$drawClasses
                            drawClasses = setNames(seq(0,1,1/((length(drawClasses)-1))), drawClasses)
                            classes <- list(
                              Inversion    = list(mult=drawClasses["Inversion"], text="HH/TT", color=self$col_h2h),
                              `Dup/Del`    = list(mult=drawClasses["Dup/Del"], text=c("DEL","DUP"), color=self$col_del),
                              Translocation= list(mult=drawClasses["Translocation"], text="TRA", color=self$col_trans)
                            )

                            for (cls in rev(names(drawClasses))) {
                              info <- classes[[cls]]
                              v_npc <- tb_y0 + (info$mult - ysc[1]) / diff(ysc) * (tb_y1 - tb_y0)

                              grid.lines(
                                x = unit(c(tb_x0, tb_x1), "npc"),
                                y = unit(rep(v_npc, 2), "npc"),
                                gp = gpar(col="grey40", lty=3, lwd=0.5)
                              )

                              if (cls == "Translocation") {
                                grid.text(
                                  label = info$text,
                                  x     = unit(tb_x0, "npc")-unit(tb_x0*0.015,"npc"),
                                  y     = unit(v_npc,"npc"),
                                  just  = "right",
                                  gp    = gpar(col="grey40", cex=0.5)
                                )
                              } else if (cls == "Dup/Del") {
                                grid.text("DEL", x=unit(tb_x0,"npc")-unit(tb_x0*0.015,"npc"),
                                          y=unit(v_npc,"npc")-unit(v_npc*0.015,"npc"),
                                          just="right", gp=gpar(col=self$col_del, cex=0.5))
                                grid.text("DUP", x=unit(tb_x0,"npc")-unit(tb_x0*0.015,"npc"),
                                          y=unit(v_npc,"npc")+unit(v_npc*0.015,"npc"),
                                          just="right", gp=gpar(col=self$col_dup, cex=0.5))
                              } else if (cls == "Inversion") {
                                grid.text("TT", x=unit(tb_x0,"npc")-unit(tb_x0*0.015,"npc"),
                                          y=unit(v_npc,"npc")-unit(v_npc*0.015,"npc"),
                                          just="right", gp=gpar(col=self$col_t2t, cex=0.5))
                                grid.text("HH", x=unit(tb_x0,"npc")-unit(tb_x0*0.015,"npc"),
                                          y=unit(v_npc,"npc")+unit(v_npc*0.015,"npc"),
                                          just="right", gp=gpar(col=self$col_h2h, cex=0.5))
                                super$draw()
                              }
                            }
                          }
                        )
)