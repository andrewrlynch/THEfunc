# SeqGene ----
#' SeqGene R6 Class
#'
#' @description
#' R6 class for plotting genes with exon–intron structures in the SeqPlot
#' framework. Each gene is drawn as a backbone line with exons as boxes,
#' directional arrows along the backbone, and an optional gene label.
#'
#' @details
#' `SeqGene` is designed for displaying annotated genes from a `GRanges`
#' input containing exon ranges. Exons are grouped by gene ID (specified
#' with `geneCol`) and arranged into non-overlapping tiers within each
#' track window. Strand information (from `strandCol` or `strand()`) is
#' used to orient arrowheads along the backbone. Labels and exon boxes
#' are automatically scaled to the track coordinate system. Colors can
#' be assigned globally (`color`) or per-gene (`colorCol`).
#'
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(
#'   seqnames = "chr1",
#'   ranges   = IRanges(c(100, 200, 400), width = 50),
#'   gene_id  = c("GENE1", "GENE1", "GENE2"),
#'   strand   = c("+", "+", "-")
#' )
#' gene_plot <- SeqGene(gr, geneCol = "gene_id", strandCol = "strand")
#' # gene_plot$prep(layout_track, track_windows)
#' # gene_plot$draw()
#'
#' @export
SeqGene <- R6::R6Class("SeqGene",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing exons with metadata.
                         gr = NULL,

                         #' @field geneCol Character scalar naming the metadata column used
                         #'   to group exons by gene (default `"gene_name"`).
                         geneCol = NULL,

                         #' @field genesFilter Optional character vector of gene IDs to keep.
                         genesFilter = NULL,

                         #' @field strandCol Optional column name giving strand orientation
                         #'   (`"+"` or `"-"`) per gene. Defaults to `strand(gr)`.
                         strandCol = NULL,

                         #' @field color Default fill/line color for genes (if `colorCol` not used).
                         color = NULL,

                         #' @field colorCol Optional column name giving per-gene colors.
                         colorCol = NULL,

                         #' @field shape Shape used for exons (currently `"rect"`).
                         shape = NULL,

                         #' @field label_pad Minimum padding (in bp) around genes to avoid
                         #'   overlapping labels. Default `50000`.
                         label_pad = 50000,

                         #' @field exon_height Proportion of tier height allocated to exon boxes.
                         exon_height = 0.8,

                         #' @field arrow_spacing Spacing between directional arrows (in bp).
                         arrow_spacing = 10,

                         #' @field arrow_len Arrow shaft length (in bp).
                         arrow_len = 1,

                         #' @field label_cex Expansion factor for gene label text. Default `0.6`.
                         label_cex = 0.6,

                         #' @field label_offset Offset for label placement (in npc units).
                         label_offset = 0.01,

                         #' @field coordCanvas Data frame of canvas coordinates built in `prep()`.
                         coordCanvas = NULL,

                         #' @description
                         #' Create a new `SeqGene` object.
                         #'
                         #' @param gr Deprecated; use `x` instead.
                         #' @param x Preferred: A `GRanges` object of exons with metadata columns for
                         #'   gene ID and (optionally) strand.
                         #' @param geneCol Metadata column giving the gene ID (default `"gene_name"`).
                         #' @param genesFilter Optional vector of gene IDs to keep.
                         #' @param strandCol Optional metadata column for strand orientation.
                         #' @param color Default color if no `colorCol` is provided.
                         #' @param colorCol Optional column for per-gene colors.
                         #' @param shape Shape for exon boxes (`"rect"` by default).
                         #' @return A new `SeqGene` object.
                         initialize = function(gr = NULL,
                                               x = NULL,
                                               geneCol     = "gene_name",
                                               genesFilter = NULL,
                                               strandCol   = NULL,
                                               color       = "gray30",
                                               colorCol    = NULL,
                                               shape       = "rect") {
                           if (!is.null(x)) {
                             stopifnot(inherits(x, "GRanges"))
                             gr_use <- x
                           } else {
                             stopifnot(inherits(gr, "GRanges"))
                             gr_use <- gr
                           }
                           self$gr          <- gr_use
                           self$geneCol     <- geneCol
                           self$genesFilter <- genesFilter
                           self$strandCol   <- strandCol
                           self$color       <- color
                           self$colorCol    <- colorCol
                           self$shape       <- shape
                           self$coordOriginal <- gr_use
                         },

                         #' @description
                         #' Prepare the canvas coordinates for plotting exons, backbones,
                         #' and labels within each track window.
                         #'
                         #' @param layout_track List of layout panels for this track.
                         #' @param track_windows A `GRanges` object of genomic windows.
                         #' @return Populates `coordCanvas` with per-exon coordinates.
                         prep = function(layout_track, track_windows) {
                           panels <- layout_track
                           nWin   <- length(panels)
                           all_rows <- list()

                           for (w in seq_len(nWin)) {
                             pm  <- panels[[w]]
                             win <- track_windows[w]

                             ws <- start(win)
                             we <- end(win)
                             if (!isTRUE(all.equal(pm$data_x, c(ws, we)))) {
                               stop(sprintf(
                                 "Window #%d mismatch: pm$data_x = [%g, %g], but window = [%g, %g]",
                                 w, pm$data_x[1], pm$data_x[2], ws, we
                               ))
                             }

                             hits <- findOverlaps(self$gr, win)
                             if (length(hits)==0) next
                             exons <- self$gr[unique(queryHits(hits))]

                             gid_all <- as.character(mcols(exons)[[self$geneCol]])
                             if (!is.null(self$genesFilter)) {
                               keep <- gid_all %in% self$genesFilter
                               exons <- exons[keep]
                               gid_all <- gid_all[keep]
                               if (length(exons)==0) next
                             }

                             mcols(exons)$gid <- gid_all
                             gene_map <- split(exons, mcols(exons)$gid)
                             gids     <- names(gene_map)

                             gene_start <- integer(length(gids))
                             gene_end   <- integer(length(gids))
                             for (i in seq_along(gids)) {
                               g <- gene_map[[gids[i]]]
                               gene_start[i] <- min(start(g))
                               gene_end[i]   <- max(end(g))
                             }

                             strand_g <- vapply(gene_map, function(g) {
                               if (!is.null(self$strandCol) && self$strandCol %in% names(mcols(g))) {
                                 s0 <- as.character(mcols(g)[[self$strandCol]][1])
                               } else {
                                 s0 <- as.character(strand(g)[1])
                               }
                               if (length(s0)!=1 || is.na(s0) || !(s0 %in% c("+","-"))) s0 <- "+"
                               s0
                             }, character(1))

                             label_g <- gids
                             if (!is.null(self$colorCol) && self$colorCol %in% names(mcols(exons))) {
                               color_g <- vapply(gene_map, function(g) mcols(g)[[self$colorCol]][1], character(1))
                             } else {
                               color_g <- setNames(rep(self$color, length(gids)), gids)
                             }

                             pad <- self$label_pad
                             pm <- panels[[w]]

                             pushViewport(viewport())
                             npc_label_widths <- as.numeric(stringWidth(label_g) * self$label_cex)
                             popViewport()

                             genomic_per_npc <- diff(pm$xscale) / (pm$inner$x1 - pm$inner$x0)
                             label_pad_bp    <- npc_label_widths * genomic_per_npc
                             pad_bp <- pmax(label_pad_bp, self$label_pad)

                             p0 <- ifelse(strand_g == "+", gene_start - pad_bp, gene_start)
                             p1 <- ifelse(strand_g == "-", gene_end   + pad_bp, gene_end)

                             ord <- order(p0)
                             ends_last <- numeric(0)
                             tiers <- integer(length(gids))
                             for (i in ord) {
                               s <- p0[i]; e <- p1[i]
                               w2 <- which(ends_last < s)
                               if (length(w2)==0) {
                                 tiers[i] <- length(ends_last)+1
                                 ends_last <- c(ends_last,e)
                               } else {
                                 tiers[i] <- w2[1]
                                 ends_last[w2[1]] <- e
                               }
                             }
                             ntiers <- max(tiers)

                             u0 <- (gene_start - pm$xscale[1]) / diff(pm$xscale)
                             u1 <- (gene_end   - pm$xscale[1]) / diff(pm$xscale)
                             x0c <- pm$inner$x0 + u0 * (pm$inner$x1 - pm$inner$x0)
                             x1c <- pm$inner$x0 + u1 * (pm$inner$x1 - pm$inner$x0)

                             track_h <- pm$inner$y1 - pm$inner$y0
                             row_h   <- track_h / ntiers
                             exon_h  <- row_h * self$exon_height
                             ymid   <- pm$inner$y0 + (tiers-1)*row_h + exon_h/2

                             for (i in seq_along(gids)) {
                               gene_id <- gids[i]
                               g       <- gene_map[[gene_id]]
                               ux0 <- (start(g) - pm$xscale[1]) / diff(pm$xscale)
                               ux1 <- (end(g)   - pm$xscale[1]) / diff(pm$xscale)
                               ex0 <- pm$inner$x0 + ux0*(pm$inner$x1 - pm$inner$x0)
                               ex1 <- pm$inner$x0 + ux1*(pm$inner$x1 - pm$inner$x0)
                               ey0 <- rep(ymid[i] - exon_h/2, length(g))
                               ey1 <- rep(ymid[i] + exon_h/2, length(g))

                               df_ex <- data.frame(
                                 gene    = gene_id,
                                 x0b     = x0c[i],
                                 x1b     = x1c[i],
                                 ymid    = ymid[i],
                                 dir     = ifelse(strand_g[i]=="-", -1, +1),
                                 label   = label_g[i],
                                 color   = color_g[i],
                                 exon_x0 = ex0,
                                 exon_x1 = ex1,
                                 exon_y0 = ey0,
                                 exon_y1 = ey1,
                                 stringsAsFactors = FALSE,
                                 row.names = NULL
                               )
                               all_rows[[ length(all_rows) + 1 ]] <- df_ex
                             }
                           }

                           self$coordCanvas <- do.call(rbind, all_rows)
                         },

                         #' @description
                         #' Draw gene backbones, exons, directional arrows, and labels
                         #' using the coordinates from `prep()`.
                         #'
                         #' @return Draws gene structures into the active grid viewport.
                         draw = function() {
                           df <- self$coordCanvas
                           if (is.null(df) || nrow(df)==0) return()

                           by_gene <- split(df, df$gene)
                           for (sub in by_gene) {
                             col  <- sub$color[1]
                             dir  <- sub$dir[1]
                             lbl  <- sub$label[1]
                             ym   <- sub$ymid[1]
                             x0b  <- sub$x0b[1]
                             x1b  <- sub$x1b[1]

                             grid.lines(
                               x = unit(c(x0b, x1b), "npc"),
                               y = unit(c(ym, ym), "npc"),
                               gp = gpar(col=col, lwd=1, lineend="butt")
                             )

                             spacing_npc <- 0.02
                             arrow_len_npc <- 0

                             x0 <- sub$x0b[1]
                             x1 <- sub$x1b[1]

                             if (dir > 0) {
                               x_start <- x0 + spacing_npc
                               x_end   <- x1 - spacing_npc
                             } else {
                               x_start <- x1 + spacing_npc
                               x_end   <- x0 - spacing_npc
                             }

                             if ((dir > 0 && x_start < x_end) || (dir < 0 && x_start > x_end)) {
                               xs <- seq(x_start, x_end, by = dir * spacing_npc)
                               if (length(xs) > 2) xs <- xs[-c(1,length(xs))]

                               for (x in xs) {
                                 grid.segments(
                                   x0 = unit(x - dir * arrow_len_npc, "npc"),
                                   x1 = unit(x, "npc"),
                                   y0 = unit(ym, "npc"),
                                   y1 = unit(ym, "npc"),
                                   gp = gpar(col = col, lwd = 1),
                                   arrow = arrow(type = "open", angle = 45, length = unit(1, "mm"))
                                 )
                               }
                             }

                             for (j in seq_len(nrow(sub))) {
                               grid.rect(
                                 x      = unit((sub$exon_x0[j]+sub$exon_x1[j])/2, "npc"),
                                 y      = unit((sub$exon_y0[j]+sub$exon_y1[j])/2, "npc"),
                                 width  = unit(sub$exon_x1[j]-sub$exon_x0[j],     "npc"),
                                 height = unit(sub$exon_y1[j]-sub$exon_y0[j],     "npc"),
                                 gp     = gpar(fill=col, col=NA)
                               )
                             }

                             labx <- if (dir>0) x0b - self$label_offset else x1b + self$label_offset
                             just <- if (dir>0) c("right","center") else c("left","center")
                             grid.text(
                               label = lbl,
                               x     = unit(labx,"npc"),
                               y     = unit(ym,"npc"),
                               just  = just,
                               gp    = gpar(cex=self$label_cex, col=col)
                             )
                           }
                         }
                       )
)