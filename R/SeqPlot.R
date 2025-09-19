#library(grid)
#library(R6)
#library(GenomicRanges)
#library(rtracklayer)

# Utilities ----

#' Load GTF file for gene tracks
#'
#' This function loads a Gencode GTF file for drawing gene tracks
#' If no path is given, it will attempt to load from the database (connection required)
#'
#' @param path Optional path to a Gencode-style gene table (GTF format).
#'
#' @return A GRange object with gene annotations
#' @export
getGencode <- function(genome = "hg38", version = "v32", path = NULL) {
  # Load from Path if given
  if(!is.null(path)){
    gtf_path <- path
    gtf <- import(gtf_path, format="gtf")
    gtf <- keepSeqlevels(gtf, paste0("chr", c(1:22,"X","Y")), pruning.mode="coarse")
    return(gtf)
  } else {
    # Load AnnotationHub
    ah <- AnnotationHub::AnnotationHub()

    # Query GENCODE GTFs matching genome and version
    query_terms <- c("GENCODE", "GTF", genome, version)
    q <- AnnotationHub::query(ah, query_terms)

    # Only keep actual GRanges GTFs, not TxDb databases
    q <- q[q$rdataclass == "GRanges"]

    if (length(q) == 0) {
      stop("No matching GENCODE GTF (GRanges) found. Check genome and version inputs.")
    }

    # Pull the first match
    gtf <- ah[[names(q)[1]]]

    # Keep only chr1-22, chrX, chrY
    gr <- GenomeInfoDb::keepSeqlevels(
      gtf,
      paste0("chr", c(1:22, "X", "Y")),
      pruning.mode = "coarse"
    )

    # Set seqlevels to UCSC style (chr-prefixed)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

    return(gr)
  }
}



#' Load cytoband table for ideogram rendering
#'
#' This function loads a UCSC-style cytoband file for ideogram drawing.
#' If no path is given, it loads the built-in cytoband_hg38 dataset from the package.
#'
#' @param path Optional path to a UCSC-style cytoband table (TSV format).
#'
#' @return A data.frame with cytoband annotations.
#' @export
loadCytobands <- function(path = NULL) {
  if (is.null(path)) {
    data("cytoband_hg38", package = "THEfunc", envir = environment())
    return(cytoband_hg38)
  } else {
    return(read.table(
      path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    ))
  }
}



#' Load colors from Kepano's Flexoki palette
#'
#' Loads colors from Kepano's Flexoki palette.
#'
#' @param n Number of colors to be produced. If > 9, an equally distributed color ramp will be produced (and probably look bad).
#'
#' @return Vector of hex colors
#' @export
flexoki_palette <- function(n, shade = 400) {
  flexoki600 <- c("#AF3029","#205EA6","#66800B","#BC5215","#24837B","#5E409D","#AD8301","#A02C6D","#6F6E69")
  flexoki400 <- c("#D14D41","#4385BE","#879A39","#DA702C","#3AA99F","#8B7EC8","#D0A215","#CE5D97","#9F9D96")
  flexoki300 <- c("#E8705F","#66A0C8","#A0AF54","#EC8B49","#5ABDAC","#A699D0","#DFB431","#E47DA8","#B7B5AC")
  flexoki200 <- c("#F89A8A","#92BFDB","#BEC97E","#F9AE77","#87D3C3","#C4B9E0","#ECCB60","#F4A4C2","#CECDC3")

  base <- if(shade == 600){flexoki600} else if (shade == 400){flexoki400} else if (shade == 300){flexoki300} else if (shade == 200){flexoki200} else {stop("shade must be 200, 300, 400, or 600")}

  if (n <= length(base)) {
    return(base[seq_len(n)])
  } else {
    # Interpolation in HCL space
    colorRampPalette(base, space = "Lab")(n)
  }
}



#' clipToXscale
#'
#' Clipping for ranges intersected by x axis limits
#'
#' @param x0,x1,xscale Coordinates for clipping
#'
#' @export
clipToXscale <- function(x0, x1, xscale) {
  keep <- x1 >= xscale[1] & x0 <= xscale[2]
  x0_clipped <- pmax(x0[keep], xscale[1])
  x1_clipped <- pmin(x1[keep], xscale[2])
  return(list(x0 = x0_clipped, x1 = x1_clipped, mask = keep))
}



#' convertDataToGrid
#'
#' Converts data coordinates to grid coordinates
#'
#' @param x,y,xscale,yscale Coordinates for clipping
#'
#' @export
convertDataToGrid <- function(x, y, xscale, yscale) {
  u <- (x - xscale[1]) / diff(xscale)
  v <- (y - yscale[1]) / diff(yscale)
  c(u, v)
}



#' gridToCanvas
#'
#' Converts grid coordinates to canvas coordinates
#'
#' @param u,v,panel_meta Coordinates for clipping
#'
#' @export
gridToCanvas <- function(u, v, panel_meta) {
  x <- panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0)
  y <- panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)
  c(x, y)
}



#' Unless A exists, then B
#'
#' @param a,b Tested quantities
#'
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b



#' defaultGenomeWindows
#' @name defaultGenomeWindows
#' @description
#' Generates a GRanges object of all standard hg38 chromosomes to serve as the default plotting windows.
#' @export
#' @author Andrew Lynch
defaultGenomeWindows <- function(add_chr = T) {
  chr_lengths <- c(
    "1"  = 248956422,
    "2"  = 242193529,
    "3"  = 198295559,
    "4"  = 190214555,
    "5"  = 181538259,
    "6"  = 170805979,
    "7"  = 159345973,
    "8"  = 145138636,
    "9"  = 138394717,
    "10" = 133797422,
    "11" = 135086622,
    "12" = 133275309,
    "13" = 114364328,
    "14" = 107043718,
    "15" = 101991189,
    "16" = 90338345,
    "17" = 83257441,
    "18" = 80373285,
    "19" = 58617616,
    "20" = 64444167,
    "21" = 46709983,
    "22" = 50818468,
    "X"  = 156040895,
    "Y"  = 57227415
  )

  if (add_chr == T){
  GRanges(seqnames = paste0("chr",names(chr_lengths)),
          ranges = IRanges(start = 1, end = chr_lengths))
  } else {
    GRanges(seqnames = names(chr_lengths),
            ranges = IRanges(start = 1, end = chr_lengths))
  }
}



#' CreateSequenceWindows
#' @name CreateSequenceWindows
#' @description
#' Generates a GRanges object of given sequence ranges
#' @param regions Vector of region strings such as "1:300-4000" or "chr2:20000-40000" or simply "chr3.
#' @param padding Number of base pairs around given positions to add. Automatically clipped to chromosome start positions. Default: 0
#' @param genome Genome of origin. Default: hg38
#' @param add_chr Whether to maintain "chr" convention on seqnames. Default: TRUE
#' @author Andrew Lynch
#' @export
CreateSequenceWindows <- function(regions, padding = 0, genome = "hg38", add_chr = TRUE) {
  regions <- gsub("chr|Chr|chromosome", "", as.character(regions))
  is_chr_only <- !grepl(":", regions)

  #Parse chrom:start-end regions
  parsed <- strcapture(
    pattern = "^([^:]+):([0-9,]+)-([0-9,]+)$",
    x = regions[!is_chr_only],
    proto = list(chrom = character(), start = character(), end = character())
  )

  #Add full-chromosome windows if applicable
  if (any(is_chr_only)) {
    chr_df <- data.frame(
      chrom = regions[is_chr_only],
      stringsAsFactors = FALSE
    )

    default_ranges <- defaultGenomeWindows(add_chr = F)
    chr_df$start <- start(default_ranges)[match(chr_df$chrom, as.character(seqnames(default_ranges)))]
    chr_df$end   <- end(default_ranges)[match(chr_df$chrom, as.character(seqnames(default_ranges)))]

    parsed <- rbind(parsed, chr_df)
  }

  parsed$start <- as.integer(parsed$start)
  parsed$end   <- as.integer(parsed$end)

  #Apply padding and clip
  parsed$start <- parsed$start - padding
  parsed$end   <- parsed$end + padding
  parsed$start[parsed$start < 1] <- 1

  #Apply chr prefix logic
  chrom <- if (add_chr) paste0("chr", parsed$chrom) else parsed$chrom

  #Build GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )

  # Merge overlapping windows
  gr_merged <- GenomicRanges::reduce(gr)

  # Return
  sort(gr_merged)
}



#' drawSeqArch
#'
#' General function for drawing SeqLink arches
#'
#' @param x0,x1,y0,y1,top0,top1 Canvas coordinates. y0 and y1 are bottoms of arch stems. top0 and top1 are the tops of arch stems
#' @param orientation Orientation of the link. In general, "-" arches will point down, "+" will point up.
#' @param curve Argument to set curve type or degree. Can be one of  "length", "equal", or number to control height of arc (default is 0.2).
#' @param stemWidth,arcWidth,arcColor,stemColor Aesthetic control over arches.
#'
#' @export
drawSeqArch <- function(x0, y0, x1, y1, top0, top1,
                        orientation = "*", curve = "length",
                        stemWidth = 1, arcWidth = 1,
                        arcColor = "black", stemColor = "black") {

  #upward <- orientation %in% c("*", "up", "+", "-/+", "-/-")
  #sign <- ifelse(upward, 1, -1)

  span <- abs(x1 - x0)

  # Determine vertical arc "bulge" height
  if (is.numeric(curve)) {
    curve_offset <- curve
  } else if (curve == "equal") {
    curve_offset <- 0.2  # default uniform curve
  } else if (curve == "length") {
    curve_offset <- span * 0.2  # scale with span (same multiplier)
  } else {
    warning("Unknown curve value in drawSeqArch(); using default 0.2")
    curve_offset <- span * 0.2
  }

  curve_offset1 <- ifelse(grepl("^\\*|^\\+", orientation), curve_offset, curve_offset * -1)
  curve_offset2 <- ifelse(grepl("\\*$|\\+$", orientation), curve_offset, curve_offset * -1)


  # Controls how broad the arch appears (smaller = broader)
  ctrl_spread <- 0

  # Horizontal offset from endpoints
  dx <- abs(x1 - x0)
  ctrl_dx <- dx * ctrl_spread

  mid_y <- max(top0, top1)
  ctrl_x1 <- x0 + ctrl_dx
  ctrl_x2 <- x1 - ctrl_dx

  ctrl_y1 <- mid_y + curve_offset1
  ctrl_y2 <- mid_y + curve_offset2

  P0 <- c(x0, top0)
  P1 <- c(ctrl_x1, ctrl_y1)
  P2 <- c(ctrl_x2, ctrl_y2)
  P3 <- c(x1, top1)

  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^3 * P0[1] + 3 * (1 - t)^2 * t * P1[1] +
    3 * (1 - t) * t^2 * P2[1] + t^3 * P3[1]
  bez_y <- (1 - t)^3 * P0[2] + 3 * (1 - t)^2 * t * P1[2] +
    3 * (1 - t) * t^2 * P2[2] + t^3 * P3[2]

  grid.lines(bez_x, bez_y, gp = gpar(col = arcColor, lwd = arcWidth))

  # Draw vertical stems
  grid.segments(x0, y0, x0, top0, gp = gpar(col = stemColor, lwd = stemWidth))
  grid.segments(x1, y1, x1, top1, gp = gpar(col = stemColor, lwd = stemWidth))
}






#' SeqElement R6 Class
#'
#' @description
#' Base R6 class representing a generic sequence element in the SeqPlot
#' framework. Stores the original genomic ranges and provides slots for
#' transformed coordinates in grid space.
#'
#' @details
#' The `SeqElement` class is the root R6 class for all genomic elements
#' in the plotting system. It is initialized with a `GRanges` object and
#' retains both the original genomic coordinates (`coordOriginal`) and
#' transformed grid coordinates (`coordGrid`).
#'
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:5, width = 1))
#' el <- SeqElement$new(gr)
#' el$gr
#'
#' @export
SeqElement <- R6::R6Class("SeqElement",
                          public = list(

                            #' @field gr A `GRanges` object containing the genomic coordinates
                            #'   and any associated metadata.
                            gr = NULL,

                            #' @field coordOriginal A `GRanges` object storing the unmodified input
                            #'   coordinates.
                            coordOriginal = NULL,

                            #' @field coordGrid Placeholder for transformed coordinates after applying
                            #'   the SeqPlot grid layout.
                            coordGrid = NULL,

                            #' @description
                            #' Create a new `SeqElement` object.
                            #' @param gr A `GRanges` object containing genomic positions and optional
                            #'   metadata columns.
                            #' @param yCol Optional character string naming a metadata column in `gr`
                            #'   to be used for y-axis values.
                            #' @return A new `SeqElement` object.
                            initialize = function(gr, yCol = NULL) {
                              if (!inherits(gr, "GRanges")) {
                                stop("gr must be a GRanges object.")
                              }
                              self$gr <- gr
                              self$coordOriginal <- gr
                            }
                          )
)




# SeqPoint ----
#' SeqPoint R6 Class
#'
#' @description
#' An R6 class for plotting genomic points, such as SNPs or single-base
#' features, on a SeqPlot track. Each genomic range is drawn as a point,
#' with optional y-axis values from metadata.
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

                          #' @field coordOriginal A `GRanges` object storing the original,
                          #'   unmodified genomic coordinates.
                          coordOriginal = NULL,

                          #' @field coordCanvas A list of per-panel coordinate matrices
                          #'   (x, y in canvas units) produced by `prep()`.
                          coordCanvas = NULL,

                          #' @field aesthetics A list of plotting aesthetics for the points,
                          #'   merged from user input and defaults.
                          aesthetics = NULL,

                          #' @field defaultAesthetics Default aesthetics for points:
                          #'   - `shape`: plotting symbol (default 16)
                          #'   - `size`: point size scaling factor
                          #'   - `color`: point color
                          defaultAesthetics = list(
                            shape = 16,
                            size = 0.1,
                            color = "#1C1B1A"
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
                          #' pt <- SeqPoint$new(gr)
                          #' pt$prep(layout_track, track_windows)
                          #' pt$draw()
                          initialize = function(gr, yCol = NULL, aesthetics = list()) {
                            stopifnot(inherits(gr, "GRanges"))
                            self$gr <- gr
                            self$coordOriginal <- gr
                            self$yCol <- yCol

                            if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
                              self$y <- as.numeric(mcols(gr)[[yCol]])
                            } else {
                              self$y <- rep(0.5, length(gr))
                            }

                            self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)
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
                              if (is.null(coords)) next
                              grid.points(
                                x = unit(coords[, 1], "npc"),
                                y = unit(coords[, 2], "npc"),
                                pch = self$aesthetics$shape,
                                gp = gpar(col = self$aesthetics$color, cex = self$aesthetics$size)
                              )
                            }
                          }
                        )
)





# SeqSegment ----
#' SeqSegment R6 Class
#'
#' @description
#' R6 class for plotting line segments in the SeqPlot framework.
#' Inherits from [SeqElement].
#'
#' @details
#' The `SeqSegment` class represents genomic features drawn as horizontal
#' or vertical line segments spanning ranges. It supports y-values from a
#' single column (`yCol`) or separate start/end columns (`y0Col`, `y1Col`).
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100), width = 50),
#'   score = c(0.2, 0.8)
#' )
#' seg <- SeqSegment$new(gr, yCol = "score")
#' seg$prep(layout_track = some_layout, track_windows = some_windows)
#' seg$draw()
#'
#' @export
SeqSegment <- R6::R6Class("SeqSegment",
                          inherit = SeqElement,
                          public = list(

                            #' @field gr A `GRanges` object containing genomic intervals.
                            gr = NULL,

                            #' @field y0 Numeric vector of lower y-values for segments.
                            y0 = NULL,

                            #' @field y1 Numeric vector of upper y-values for segments.
                            y1 = NULL,

                            #' @field yCol Optional metadata column used for both y0 and y1 values.
                            yCol = NULL,

                            #' @field y0Col Optional metadata column used for lower y-values.
                            y0Col = NULL,

                            #' @field y1Col Optional metadata column used for upper y-values.
                            y1Col = NULL,

                            #' @field coordOriginal A `GRanges` object storing the unmodified input
                            #'   coordinates.
                            coordOriginal = NULL,

                            #' @field coordCanvas A list of transformed segment coordinates
                            #'   prepared for plotting.
                            coordCanvas = NULL,

                            #' @field aesthetics List of aesthetics merged with defaults.
                            aesthetics = NULL,

                            #' @field defaultAesthetics Default aesthetics:
                            #'   \code{lwd = 1.5}, \code{col = "#1C1B1A"}.
                            defaultAesthetics = list(
                              lwd = 1.5,
                              col = "#1C1B1A"
                            ),

                            #' @description
                            #' Create a new `SeqSegment` object.
                            #' @param gr A `GRanges` object containing genomic intervals.
                            #' @param yCol Optional column name in `gr` used for both y0 and y1.
                            #' @param y0Col Optional column name in `gr` for lower y-values.
                            #' @param y1Col Optional column name in `gr` for upper y-values.
                            #' @param aesthetics Optional list of aesthetic overrides.
                            #' @return A new `SeqSegment` object.
                            initialize = function(gr, yCol = NULL, y0Col = NULL, y1Col = NULL,
                                                  aesthetics = list()) {
                              stopifnot(inherits(gr, "GRanges"))
                              self$gr <- gr
                              self$coordOriginal <- gr
                              self$yCol <- yCol
                              self$y0Col <- y0Col
                              self$y1Col <- y1Col

                              if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
                                self$y0 <- self$y1 <- as.numeric(mcols(gr)[[yCol]])
                              } else {
                                if (!is.null(y0Col) && y0Col %in% names(mcols(gr))) {
                                  self$y0 <- as.numeric(mcols(gr)[[y0Col]])
                                } else {
                                  self$y0 <- rep(0.5, length(gr))
                                }

                                if (!is.null(y1Col) && y1Col %in% names(mcols(gr))) {
                                  self$y1 <- as.numeric(mcols(gr)[[y1Col]])
                                } else {
                                  self$y1 <- self$y0
                                }
                              }

                              self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)
                            },

                            #' @description
                            #' Prepare segment coordinates by mapping genomic intervals into
                            #' panel-relative canvas space.
                            #' @param layout_track A list of panel layout metadata.
                            #' @param track_windows A `GRanges` object of genomic windows.
                            prep = function(layout_track, track_windows) {
                              self$coordCanvas <- vector("list", length(track_windows))
                              ...
                            },

                            #' @description
                            #' Draw line segments on the plotting canvas.
                            draw = function() {
                              ...
                            }
                          )
)





# SeqRect ----
#' SeqRect R6 Class
#'
#' @description
#' R6 class for drawing rectangular genomic features (boxes) in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' `SeqRect` is useful for visualizing genomic intervals as filled rectangles,
#' for example in bar plots, ideograms, or feature blocks. By default, each
#' rectangle is drawn centered on a y-value (`yCol` or a constant) with a
#' configurable relative height (`width` aesthetic).
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100), width = 50),
#'   score = c(0.2, 0.8)
#' )
#' rects <- SeqRect$new(gr, yCol = "score")
#' rects$prep(layout_track = some_layout, track_windows = some_windows)
#' rects$draw()
#'
#' @export
SeqRect <- R6::R6Class("SeqRect",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic intervals.
                         gr = NULL,

                         #' @field y0 Numeric vector of lower rectangle y-values (computed in prep).
                         y0 = NULL,

                         #' @field y1 Numeric vector of upper rectangle y-values (computed in prep).
                         y1 = NULL,

                         #' @field y Numeric vector of rectangle center y-values.
                         y = NULL,

                         #' @field yCol Optional column name in `gr` used for y-values.
                         yCol = NULL,

                         #' @field coordCanvas A list of matrices containing transformed rectangle
                         #'   coordinates for each window (x0, x1, y0, y1).
                         coordCanvas = NULL,

                         #' @field aesthetics List of current aesthetics merged with defaults.
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default aesthetics: \code{fill = "grey80"},
                         #'   \code{col = "#1C1B1A"}, \code{lwd = 0.5}, \code{width = 0.1}.
                         defaultAesthetics = list(
                           fill = "grey80",
                           col = "#1C1B1A",
                           lwd = 0.5,
                           width = 0.1
                         ),

                         #' @description
                         #' Create a new `SeqRect` object.
                         #'
                         #' @param gr A `GRanges` object containing genomic intervals.
                         #' @param yCol Optional column name in `gr` used for rectangle y-centers.
                         #' @param aesthetics Optional list of aesthetic overrides.
                         #' @return A new `SeqRect` object.
                         initialize = function(gr, yCol = NULL, aesthetics = list()) {
                           stopifnot(inherits(gr, "GRanges"))
                           self$gr <- gr
                           self$yCol <- yCol
                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)

                           y_center <- if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
                             as.numeric(mcols(gr)[[yCol]])
                           } else {
                             rep(0.5, length(gr))
                           }

                           self$y <- y_center
                         },

                         #' @description
                         #' Prepare rectangle coordinates by mapping genomic intervals into
                         #' panel-relative canvas space.
                         #'
                         #' @param layout_track A list of panel layout metadata.
                         #' @param track_windows A `GRanges` object of genomic windows.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))
                           ov <- findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- queryHits(ov)
                           sh <- subjectHits(ov)
                           x0 <- start(self$gr)[qh]
                           x1 <- end(self$gr)[qh]
                           y_center <- self$y[qh]

                           for (w in unique(sh)) {
                             mask <- sh == w
                             p <- layout_track[[w]]

                             # Subset to current panel's hits
                             x0_sub <- x0[mask]
                             x1_sub <- x1[mask]
                             y_sub <- y_center[mask]

                             # Clip horizontally
                             clip <- clipToXscale(x0_sub, x1_sub, p$xscale)

                             if (length(clip$x0) == 0) {
                               self$coordCanvas[[w]] <- NULL
                               next
                             }

                             x0_sub <- clip$x0
                             x1_sub <- clip$x1
                             y_sub  <- y_sub[clip$mask]

                             # Sanity check
                             if (length(y_sub) == 0) {
                               self$coordCanvas[[w]] <- matrix(
                                 numeric(0),
                                 ncol = 4,
                                 dimnames = list(NULL, c("x0", "x1", "y0", "y1"))
                               )
                               next
                             }

                             # Transform X to panel-relative [0–1]
                             u0 <- (x0_sub - p$xscale[1]) / diff(p$xscale)
                             u1 <- (x1_sub - p$xscale[1]) / diff(p$xscale)

                             # Vertical position and height
                             v_center <- (y_sub - p$yscale[1]) / diff(p$yscale)
                             v_center <- pmax(pmin(v_center, 1), 0)
                             v_half <- self$aesthetics$width / 2
                             v0 <- v_center - v_half
                             v1 <- v_center + v_half

                             # Convert to canvas space
                             xleft <- p$inner$x0 + u0 * (p$inner$x1 - p$inner$x0)
                             xright <- p$inner$x0 + u1 * (p$inner$x1 - p$inner$x0)
                             ybottom <- p$inner$y0 + v0 * (p$inner$y1 - p$inner$y0)
                             ytop <- p$inner$y0 + v1 * (p$inner$y1 - p$inner$y0)

                             # Store as matrix for vectorized drawing
                             self$coordCanvas[[w]] <- cbind(
                               x0 = xleft,
                               x1 = xright,
                               y0 = ybottom,
                               y1 = ytop
                             )
                           }
                         },

                         #' @description
                         #' Draw rectangles on the plotting canvas.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()

                           for (coords in self$coordCanvas) {
                             if (!is.matrix(coords) || nrow(coords) == 0 || is.null(colnames(coords))) next

                             grid.rect(
                               x = unit((coords[, "x0"] + coords[, "x1"]) / 2, "npc"),
                               y = unit((coords[, "y0"] + coords[, "y1"]) / 2, "npc"),
                               width = unit(abs(coords[, "x1"] - coords[, "x0"]), "npc"),
                               height = unit(abs(coords[, "y1"] - coords[, "y0"]), "npc"),
                               gp = gpar(
                                 fill = self$aesthetics$fill,
                                 col = self$aesthetics$col,
                                 lwd = self$aesthetics$lwd
                               )
                             )
                           }
                         }
                       )
)





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
#' bars <- SeqBar$new(gr, yCol = "score", groupCol = "group")
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
                        initialize = function(gr, yCol = NULL, groupCol = NULL,
                                              groupLevels = NULL, aesthetics = list()) {
                          stopifnot(inherits(gr, "GRanges"))
                          self$gr <- gr
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





# SeqLine ----
#' SeqLine R6 Class
#'
#' @description
#' R6 class for drawing line plots from genomic intervals in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' Each genomic interval is represented as a point at its midpoint by default,
#' with y-values taken from a metadata column or set to a constant. The points
#' are connected into a continuous line. Step lines can also be drawn by setting
#' the `type` aesthetic to `"s"` or `"step"`.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   "chr1",
#'   IRanges::IRanges(c(1, 100, 200), width = 50),
#'   score = c(2, 5, 3)
#' )
#' line <- SeqLine$new(gr, yCol = "score")
#' line$prep(layout_track = some_layout, track_windows = some_windows)
#' line$draw()
#'
#' @export
SeqLine <- R6::R6Class("SeqLine",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr A `GRanges` object containing genomic intervals.
                         gr = NULL,

                         #' @field y Numeric vector of y-values (from `yCol` or constant).
                         y = NULL,

                         #' @field yCol Optional column name in `gr` used for y-values.
                         yCol = NULL,

                         #' @field coordOriginal A `GRanges` object storing unmodified input coordinates.
                         coordOriginal = NULL,

                         #' @field coordCanvas List of matrices storing transformed line coordinates
                         #'   in canvas space for each genomic window.
                         coordCanvas = NULL,

                         #' @field aesthetics List of current aesthetics merged with defaults.
                         aesthetics = NULL,

                         #' @field defaultAesthetics Default aesthetics for lines:
                         #'   \code{type = "n"}, \code{size = 0.1}, \code{color = "#1C1B1A"}.
                         defaultAesthetics = list(
                           type = "n",
                           size = 0.1,
                           color = "#1C1B1A"
                         ),

                         #' @description
                         #' Create a new `SeqLine` object.
                         #'
                         #' @param gr A `GRanges` object containing genomic intervals.
                         #' @param yCol Optional column name in `gr` for y-values.
                         #' @param aesthetics Optional list of aesthetic overrides.
                         #' @return A new `SeqLine` object.
                         initialize = function(gr, yCol = NULL, aesthetics = list()) {
                           stopifnot(inherits(gr, "GRanges"))
                           self$gr <- gr
                           self$coordOriginal <- gr
                           self$yCol <- yCol

                           if (!is.null(yCol) && yCol %in% names(mcols(gr))) {
                             self$y <- as.numeric(mcols(gr)[[yCol]])
                           } else {
                             self$y <- rep(0.5, length(gr))
                           }

                           self$aesthetics <- modifyList(self$defaultAesthetics, aesthetics)
                         },

                         #' @description
                         #' Prepare line coordinates in canvas space for each genomic window.
                         #'
                         #' @param layout_track A list of panel layout metadata for a track.
                         #' @param track_windows A `GRanges` object defining genomic windows.
                         prep = function(layout_track, track_windows) {
                           self$coordCanvas <- vector("list", length(track_windows))

                           ov <- findOverlaps(self$gr, track_windows)
                           if (length(ov) == 0) return(invisible())

                           qh <- S4Vectors::queryHits(ov)
                           sh <- S4Vectors::subjectHits(ov)

                           x <- (start(self$gr)[qh] + end(self$gr)[qh]) / 2
                           y <- self$y[qh]

                           if (self$aesthetics$type %in% c("s", "step")) {
                             x <- start(self$gr)[qh]
                             x <- rep(x, each = 2)[-1]
                             y <- rep(y, each = 2)[-length(y) * 2]
                           }

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
                           }

                           invisible()
                         },

                         #' @description
                         #' Draw lines onto the plotting canvas.
                         draw = function() {
                           if (is.null(self$coordCanvas)) return()
                           for (w in seq_along(self$coordCanvas)) {
                             coords <- self$coordCanvas[[w]]
                             if (is.null(coords)) next
                             grid.lines(
                               x = unit(coords[, 1], "npc"),
                               y = unit(coords[, 2], "npc"),
                               gp = gpar(col = self$aesthetics$color, cex = self$aesthetics$size)
                             )
                           }
                         }
                       )
)





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
#' area <- SeqArea$new(gr, yCol = "score", groupCol = "group")
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
                         initialize = function(gr, yCol = NULL, groupCol = NULL,
                                               groupLevels = NULL, aesthetics = list()) {
                           stopifnot(inherits(gr, "GRanges"))
                           self$gr <- gr
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


# SeqLink ----
#' SeqLink R6 Class
#'
#' @description
#' R6 class representing a link between two genomic loci in the SeqPlot
#' framework. Inherits from [SeqElement].
#'
#' @details
#' A `SeqLink` connects two sets of genomic intervals (`gr1` and `gr2`),
#' optionally across different tracks. It stores track indices (`t0`, `t1`),
#' vertical anchor positions (`y0`, `y1`), and link orientations. Higher-level
#' link classes such as `SeqArch` extend this class to implement specific
#' visualization styles (arches, lines, bezier curves, etc.).
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 50), width = 1))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), width = 1))
#' link <- SeqLink$new(gr1, gr2, t0 = 1, t1 = 1, y0 = 0, y1 = 0, orientation = "+")
#' link$gr1
#'
#' @export
SeqLink <- R6::R6Class("SeqLink",
                       inherit = SeqElement,
                       public = list(

                         #' @field gr1 A `GRanges` object containing the start loci of each link.
                         gr1 = NULL,

                         #' @field gr2 A `GRanges` object containing the end loci of each link.
                         gr2 = NULL,

                         #' @field t0 Integer vector giving the track indices of the start loci.
                         t0 = NULL,

                         #' @field t1 Integer vector giving the track indices of the end loci.
                         t1 = NULL,

                         #' @field y0 Numeric vector of vertical anchor positions for the start loci.
                         y0 = NULL,

                         #' @field y1 Numeric vector of vertical anchor positions for the end loci.
                         y1 = NULL,

                         #' @field orientation Character vector specifying orientation of each link
                         #'   (e.g. `"+"`, `"-"`, `"*"`).
                         orientation = NULL,

                         #' @description
                         #' Create a new `SeqLink` object.
                         #'
                         #' @param gr1 A `GRanges` object giving start positions of links.
                         #' @param gr2 A `GRanges` object giving end positions of links.
                         #' @param t0 Integer or vector giving track indices for start loci (default: `0`).
                         #' @param t1 Integer or vector giving track indices for end loci (default: `0`).
                         #' @param y0 Numeric or vector giving vertical start positions (default: `0`).
                         #' @param y1 Numeric or vector giving vertical end positions (default: `0`).
                         #' @param color Color for links (currently unused in base class).
                         #' @param orientation Character or vector specifying link orientation
                         #'   (default: `"*"`).
                         #' @return A new `SeqLink` object.
                         initialize = function(gr1, gr2, t0 = 0, t1 = 0,
                                               y0 = 0, y1 = 0,
                                               color = "black", orientation = "*") {
                           stopifnot(length(gr1) == length(gr2))
                           self$gr1 <- gr1
                           self$gr2 <- gr2
                           self$t0 <- rep(t0, length(gr1))
                           self$t1 <- rep(t1, length(gr1))
                           self$y0 <- rep(y0, length(gr1))
                           self$y1 <- rep(y1, length(gr2))
                           self$orientation <- rep(orientation, length(gr1))
                         }
                       )
)




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
#' arch <- SeqArch$new(gr1, gr2, t0 = 1, t1 = 1, height = 0.5)
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
                         #' @param gr1 A `GRanges` object giving start positions of links.
                         #' @param gr2 A `GRanges` object giving end positions of links.
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
                         initialize = function(gr1, gr2, t0 = 0, t1 = 0,
                                               y0 = 0, y1 = 0, yCol = NULL, orientation = "*",
                                               height = 1, curve = "length", aesthetics = list()) {

                           super$initialize(gr1 = gr1, gr2 = gr2, t0 = t0, t1 = t1, y0 = y0, y1 = y1,
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



# SeqIdeogram ----
#' SeqIdeogram R6 Class
#'
#' @description
#' An R6 class for drawing chromosome ideograms from cytogenetic band
#' annotations. Each band is drawn as a rectangle, and centromeres
#' (`acen` regions) are rendered as paired red triangles.
#'
#' @export
SeqIdeogram <- R6::R6Class("SeqIdeogram",
                           inherit = SeqElement,
                           public = list(

                             #' @field cytobands A `GRanges` object of cytogenetic bands. Must
                             #'   contain a `gieStain` metadata column.
                             cytobands = NULL,

                             #' @field coordCanvas A list of data frames holding canvas
                             #'   coordinates for non-centromeric bands (`x0`, `x1`, `y0`, `y1`)
                             #'   and their fill colors.
                             coordCanvas = NULL,

                             #' @field centroPolys A list of polygons representing centromere
                             #'   regions, each stored as x/y coordinate vectors.
                             centroPolys = NULL,

                             #' @description
                             #' Create a new `SeqIdeogram` object.
                             #' @param cytobands A `GRanges` object of cytobands with a `gieStain`
                             #'   metadata column.
                             #' @return A new `SeqIdeogram` object.
                             #' @examples
                             #' cb <- GenomicRanges::GRanges(
                             #'   seqnames = Rle("chr1"),
                             #'   ranges   = IRanges::IRanges(c(1, 500001, 1000001), width = 500000),
                             #'   gieStain = c("gneg", "acen", "gpos100")
                             #' )
                             #' ideogram <- SeqIdeogram$new(cb)
                             initialize = function(cytobands) {
                               stopifnot(inherits(cytobands, "GRanges"))
                               self$cytobands <- cytobands
                             },

                             #' @description
                             #' Prepare ideogram band and centromere coordinates for each panel
                             #' window in the SeqPlot layout.
                             #' @param layout_track A list of panel metadata from
                             #'   `SeqPlot$layoutGrid()`.
                             #' @param track_windows A `GRanges` object of genomic windows for the
                             #'   track.
                             prep = function(layout_track, track_windows) {
                               n_panels <- length(track_windows)
                               self$coordCanvas <- vector("list", n_panels)
                               self$centroPolys <- vector("list", n_panels)

                               ov <- findOverlaps(self$cytobands, track_windows)
                               if (length(ov) == 0) return()

                               qh <- queryHits(ov); sh <- subjectHits(ov)

                               for (w in unique(sh)) {
                                 panel <- layout_track[[w]]
                                 bands <- self$cytobands[qh[sh == w]]

                                 u0 <- (start(bands) - panel$xscale[1]) / diff(panel$xscale)
                                 u1 <- (end(bands)   - panel$xscale[1]) / diff(panel$xscale)
                                 u0 <- pmax(pmin(u0, 1), 0); u1 <- pmax(pmin(u1, 1), 0)

                                 x0c <- panel$inner$x0 + u0 * (panel$inner$x1 - panel$inner$x0)
                                 x1c <- panel$inner$x0 + u1 * (panel$inner$x1 - panel$inner$x0)
                                 y0c <- panel$inner$y0
                                 y1c <- panel$inner$y1

                                 stain <- mcols(bands)$gieStain
                                 fillCols <- sapply(stain, function(s) {
                                   if (s == "gneg") "#FFFFFF"
                                   else if (startsWith(s, "gpos")) {
                                     pct <- as.numeric(sub("gpos", "", s)) / 100
                                     grey(pct)
                                   } else if (s == "acen") "#FF0000"
                                   else "#CCCCCC"
                                 })

                                 non_acen <- stain != "acen"
                                 self$coordCanvas[[w]] <- data.frame(
                                   x0 = x0c[non_acen], x1 = x1c[non_acen],
                                   y0 = y0c,          y1 = y1c,
                                   fill = fillCols[non_acen],
                                   stringsAsFactors = FALSE
                                 )

                                 acen_idx <- which(stain == "acen")
                                 if (length(acen_idx) == 2) {
                                   x0a <- x0c[acen_idx[1]]; x1a <- x1c[acen_idx[1]]
                                   x0b <- x0c[acen_idx[2]]; x1b <- x1c[acen_idx[2]]
                                   ym  <- (y0c + y1c) / 2

                                   tri1_x <- c(x0a, x1a, x0a)
                                   tri1_y <- c(y0c, ym,  y1c)

                                   tri2_x <- c(x1b, x0b, x1b)
                                   tri2_y <- c(y0c, ym,  y1c)

                                   self$centroPolys[[w]] <- list(
                                     list(x = tri1_x, y = tri1_y),
                                     list(x = tri2_x, y = tri2_y)
                                   )
                                 }
                               }
                             },

                             #' @description
                             #' Draw the ideogram for all windows: non-centromeric bands as
                             #' rectangles, centromeres as paired red triangles.
                             draw = function() {
                               if (!is.null(self$coordCanvas)) {
                                 for (coords in self$coordCanvas) {
                                   if (nrow(coords) == 0) next
                                   grid.rect(
                                     x      = unit((coords$x0 + coords$x1) / 2, "npc"),
                                     y      = unit((coords$y0 + coords$y1) / 2, "npc"),
                                     width  = unit(coords$x1 - coords$x0, "npc"),
                                     height = unit(coords$y1 - coords$y0, "npc"),
                                     gp     = gpar(fill = coords$fill, col = "black", lwd = 0.66)
                                   )
                                 }
                               }
                               if (!is.null(self$centroPolys)) {
                                 for (polys in self$centroPolys) {
                                   if (is.null(polys)) next
                                   for (tri in polys) {
                                     grid.polygon(
                                       x = unit(tri$x, "npc"),
                                       y = unit(tri$y, "npc"),
                                       gp = gpar(fill = "#FF0000", col = "black", lwd = 0.66)
                                     )
                                   }
                                 }
                               }
                             }
                           )
)


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
#' recon <- SeqRecon$new(gr1, gr2)
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
                          #' @param gr1 A `GRanges` object with start loci.
                          #' @param gr2 A `GRanges` object with end loci.
                          #' @param yCol Optional column name for custom arch heights.
                          #' @param y0,y1 Numeric start and end y-values (default: `0`).
                          #' @param t0,t1 Track indices for start and end loci (default: `0` → current track).
                          #' @param orientation Character vector of orientations (default: `"*"`).
                          #' @param curve Character curvature mode (`"length"`, `"equal"`, or numeric).
                          #' @param drawClasses Which class tiers to render.
                          #' @param aesthetics List of aesthetic overrides (colors, widths, etc.).
                          #' @return A new `SeqRecon` object.
                          initialize = function(gr1, gr2,
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



# SeqRecon <- R6::R6Class("SeqRecon",
#                         inherit = SeqArch,
#                         public = list(
#                           last_arc_track = NULL,
#
#                           # new color fields
#                           col_h2h = NULL,
#                           col_t2t = NULL,
#                           col_dup = NULL,
#                           col_del = NULL,
#                           col_trans = NULL,
#                           drawClasses = c("Inversion", "Dup/Del", "Translocation"),
#
#                           arc_track_idx     = NULL,
#                           layout_all_tracks = NULL,
#
#                           tierMultipliers = NULL,
#
#                           initialize = function(gr1, gr2,
#                                                 yCol = NULL,
#                                                 y0 = 0, y1 = 0,
#                                                 t0 = 0, t1 = 0,
#                                                 orientation = "*",
#                                                 curve = "length",
#                                                 drawClasses = c("Inversion", "Dup/Del", "Translocation"),
#                                                 aesthetics = list()) {
#
#                             self$col_h2h   <- aesthetics$h2hColor   %||% flexoki_palette(9)[3]
#                             self$col_t2t   <- aesthetics$t2tColor   %||% flexoki_palette(9)[4]
#                             self$col_dup   <- aesthetics$dupColor   %||% flexoki_palette(9)[1]
#                             self$col_del   <- aesthetics$delColor   %||% flexoki_palette(9)[2]
#                             self$col_trans <- aesthetics$transColor %||% flexoki_palette(9)[9]
#
#                             super$initialize(
#                               gr1        = gr1,
#                               gr2        = gr2,
#                               t0         = t0,
#                               t1         = t1,
#                               y0         = 0,
#                               y1         = 0,
#                               yCol       = yCol,
#                               orientation= orientation,
#                               height     = 1,
#                               curve      = curve,
#                               aesthetics = aesthetics
#                             )
#
#                             self$aesthetics <- modifyList(self$aesthetics, list(drawClasses = self$drawClasses))
#                           },
#
#                           prep = function(layout_all_tracks, track_windows_list, arc_track_idx) {
#                             self$last_arc_track <- layout_all_tracks[[arc_track_idx]]
#
#                             self$coordGrid <- NULL
#
#                             N    <- length(self$gr1)
#                             seq1 <- as.character(seqnames(self$gr1))
#                             seq2 <- as.character(seqnames(self$gr2))
#
#                             s1   <- ifelse(as.character(strand(self$gr1)) %in% c("+","-"),
#                                            as.character(strand(self$gr1)), "+")
#                             s2   <- ifelse(as.character(strand(self$gr2)) %in% c("+","-"),
#                                            as.character(strand(self$gr2)), "+")
#                             code <- paste0(s1, "/", s2)
#
#                             tier  <- numeric(N)
#                             cols  <- character(N)
#                             ori   <- character(N)
#
#                             t0 <- ifelse(self$t0 == 0, arc_track_idx, self$t0)
#                             t1 <- ifelse(self$t1 == 0, arc_track_idx, self$t1)
#
#                             for (i in seq_len(N)) {
#                               if (seq1[i] != seq2[i]) {
#                                 tier[i] <- 1
#                                 cols[i] <- self$col_trans
#                                 ori[i]  <- code[i]
#                               } else if (code[i] %in% c("+/+", "-/-")) {
#                                 tier[i] <- 0
#                                 cols[i] <- if (code[i] == "+/+") self$col_h2h else self$col_t2t
#                                 ori[i]  <- ifelse(code[i] == "+/+", "+", "-")
#                               } else if (code[i] %in% c("-/+", "+/-")) {
#                                 tier[i] <- 0.5
#                                 cols[i] <- if (code[i] == "-/+") self$col_dup else self$col_del
#                                 ori[i]  <- ifelse(code[i] == "-/+", "+", "-")
#                               } else {
#                                 tier[i] <- 1
#                                 cols[i] <- self$col_trans
#                                 ori[i]  <- "+"
#                               }
#                             }
#
#                             self$tierMultipliers <- tier
#                             self$orientation <- ori
#                             self$height      <- tier
#                             self$aesthetics$arcColor  <- cols
#                             self$aesthetics$stemColor <- cols
#                             self$layout_all_tracks <- layout_all_tracks
#                             self$arc_track_idx     <- arc_track_idx
#
#                             super$prep(layout_all_tracks, track_windows_list, arc_track_idx)
#                           },
#
#                           draw = function() {
#                             # nothing to do if prep never ran
#                             if (is.null(self$last_arc_track)) return()
#
#                             # 1) compute the NPC bbox of the entire recon‐track
#                             panels <- self$last_arc_track
#                             x0s <- vapply(panels, function(pm) pm$inner$x0, numeric(1))
#                             x1s <- vapply(panels, function(pm) pm$inner$x1, numeric(1))
#                             y0s <- vapply(panels, function(pm) pm$inner$y0, numeric(1))
#                             y1s <- vapply(panels, function(pm) pm$inner$y1, numeric(1))
#                             tb_x0 <- min(x0s); tb_x1 <- max(x1s)
#                             tb_y0 <- min(y0s); tb_y1 <- max(y1s)
#                             ysc   <- panels[[1]]$yscale
#
#                             # Class tiers
#                             drawClasses = self$drawClasses
#                             drawClasses = setNames(seq(0,1,1/((length(drawClasses)-1))), drawClasses)
#                             classes <- list(
#                               Inversion    = list(mult=drawClasses["Inversion"], text="HH/TT",
#                                                   color = self$col_h2h),
#                               `Dup/Del`    = list(mult=drawClasses["Dup/Del"], text=c("DEL","DUP"),
#                                                   color = self$col_del),
#                               Translocation= list(mult=drawClasses["Translocation"], text="TRA",
#                                                   color = self$col_trans)
#                             )
#
#                             # Draw tiers
#                             for (cls in rev(names(drawClasses))) {
#                               info <- classes[[cls]]
#                               v_npc <- tb_y0 + (info$mult - ysc[1]) / diff(ysc) * (tb_y1 - tb_y0)
#
#                               grid.lines(
#                                 x = unit(c(tb_x0, tb_x1), "npc"),
#                                 y = unit(rep(v_npc, 2), "npc"),
#                                 gp = gpar(col="grey40", lty=3, lwd = 0.5)
#                               )
#
#                               # Top tier
#                               if (cls == "Translocation") {
#                                 grid.text(
#                                   label = info$text,
#                                   x     = unit(tb_x0, "npc") - unit(2, "mm"),
#                                   y     = unit(v_npc, "npc"),
#                                   just  = "right",
#                                   gp    = gpar(col="grey40", cex=0.5)
#                                 )
#                               }
#
#                               # Middle tier
#                               else if (cls == "Dup/Del") {
#                                 # DEL
#                                 grid.text(
#                                   label = info$text[1],
#                                   x     = unit(tb_x0, "npc") - unit(tb_x0*0.015, "npc"),
#                                   y     = unit(v_npc, "npc") - unit(v_npc*0.015, "npc"),
#                                   just  = "right",
#                                   gp    = gpar(col = self$col_del, cex=0.5)
#                                 )
#                                 # DUP
#                                 grid.text(
#                                   label = info$text[2],
#                                   x     = unit(tb_x0, "npc") - unit(tb_x0*0.015, "npc"),
#                                   y     = unit(v_npc, "npc") + unit(v_npc*0.015, "npc"),
#                                   just  = "right",
#                                   gp    = gpar(col = self$col_dup, cex=0.5)
#                                 )
#                               }
#
#                               # Bottom tier
#                               else if (cls == "Inversion") {
#                                 # TT
#                                 grid.text(
#                                   label = "TT",
#                                   x     = unit(tb_x0, "npc") - unit(tb_x0*0.015, "npc"),
#                                   y     = unit(v_npc, "npc") - unit(v_npc*0.015, "npc"),
#                                   just  = "right",
#                                   gp    = gpar(col = self$col_t2t, cex=0.5)
#                                 )
#                                 # HH
#                                 grid.text(
#                                   label = "HH",
#                                   x     = unit(tb_x0, "npc") - unit(tb_x0*0.015, "npc"),
#                                   y     = unit(v_npc, "npc") + unit(v_npc*0.015, "npc"),
#                                   just  = "right",
#                                   gp    = gpar(col = self$col_h2h, cex=0.5)
#                                 )
#                                 super$draw()
#
#                               }
#                             }
#
#                             # 5) lastly, draw the arches on top
#                             #super$draw()
#                           }
#
#
#
#                         )
# )



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
#' gene_plot <- SeqGene$new(gr, geneCol = "gene_id", strandCol = "strand")
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
                         #' @param gr A `GRanges` object of exons with metadata columns for
                         #'   gene ID and (optionally) strand.
                         #' @param geneCol Metadata column giving the gene ID (default `"gene_name"`).
                         #' @param genesFilter Optional vector of gene IDs to keep.
                         #' @param strandCol Optional metadata column for strand orientation.
                         #' @param color Default color if no `colorCol` is provided.
                         #' @param colorCol Optional column for per-gene colors.
                         #' @param shape Shape for exon boxes (`"rect"` by default).
                         #' @return A new `SeqGene` object.
                         initialize = function(gr,
                                               geneCol     = "gene_name",
                                               genesFilter = NULL,
                                               strandCol   = NULL,
                                               color       = "gray30",
                                               colorCol    = NULL,
                                               shape       = "rect") {
                           stopifnot(inherits(gr,"GRanges"))
                           self$gr          <- gr
                           self$geneCol     <- geneCol
                           self$genesFilter <- genesFilter
                           self$strandCol   <- strandCol
                           self$color       <- color
                           self$colorCol    <- colorCol
                           self$shape       <- shape
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

                             ws <- start(win); we <- end(win)
                             if (!identical(pm$xscale, c(ws, we))) {
                               stop(sprintf(
                                 "Window #%d mismatch: pm$xscale = [%g, %g], but window = [%g, %g]",
                                 w, pm$xscale[1], pm$xscale[2], ws, we
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





# SeqTrack ----
#' SeqTrack R6 Class
#'
#' @description
#' An R6 class representing a single track in the SeqPlot framework.
#' A track contains one or more sequence elements (`SeqElement` objects),
#' associated genomic windows (`GRanges`), and track-specific aesthetics
#' such as axis visibility and titles.
#'
#' @details
#' Tracks are the core containers in a `SeqPlot` object. Each track defines
#' its own genomic windows and holds a list of plotting elements (such as
#' `SeqPoint`, `SeqBar`, `SeqLine`, or `SeqLink`). Track-level aesthetics
#' allow customization of axes, labels, and y-axis ranges.
#'
#' @examples
#' library(GenomicRanges)
#' win <- GRanges("chr1", IRanges(c(1, 1001), width = 500))
#' track <- SeqTrack$new(windows = win)
#' track$addElement(SeqRect$new(win))
#' length(track$elements)
#'
#' @export
SeqTrack <- R6::R6Class("SeqTrack",
                        public = list(

                          #' @field elements A list of `SeqElement` objects contained in this track.
                          elements = NULL,

                          #' @field windows A `GRanges` object defining the genomic windows
                          #'   covered by this track.
                          windows = NULL,

                          #' @field aesthetics A named list of aesthetics controlling axis
                          #'   rendering and labeling:
                          #'   \describe{
                          #'     \item{xAxisTitle}{Logical; whether to show the x-axis title.}
                          #'     \item{yAxisTitle}{Logical; whether to show the y-axis title.}
                          #'     \item{yAxisTitleText}{Optional custom text for the y-axis title.}
                          #'     \item{yAxisLimits}{Optional numeric vector of length 2 setting
                          #'     the y-axis limits.}
                          #'   }
                          aesthetics = list(
                            xAxisTitle = TRUE,
                            yAxisTitle = TRUE,
                            yAxisTitleText = NULL,
                            yAxisLimits = NULL
                          ),

                          #' @description
                          #' Create a new `SeqTrack` object.
                          #'
                          #' @param elements Optional list of `SeqElement` objects to add
                          #'   when constructing the track.
                          #' @param windows Optional `GRanges` object defining the genomic
                          #'   windows for this track.
                          #' @param aesthetics Optional named list of track aesthetics,
                          #'   overriding defaults.
                          #' @return A new `SeqTrack` object.
                          initialize = function(elements = list(),
                                                windows = NULL,
                                                aesthetics = list(
                                                  xAxisTitle = TRUE,
                                                  yAxisTitle = TRUE,
                                                  yAxisTitleText = NULL,
                                                  yAxisLimits = NULL
                                                )) {
                            self$elements <- elements
                            self$windows <- windows
                            self$aesthetics <- aesthetics
                          },

                          #' @description
                          #' Add a `SeqElement` object to this track.
                          #'
                          #' @param feature A `SeqElement` object to be appended to the track.
                          #' @return Updates the `elements` field with the new feature.
                          addElement = function(feature) {
                            self$elements <- append(self$elements, list(feature))
                          }
                        )
)



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
#' track1 <- SeqTrack$new(windows = win)
#' track2 <- SeqTrack$new(windows = win)
#' sp <- SeqPlot$new(tracks = list(track1, track2), windows = win)
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
                       yAxisPerWindow = FALSE
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
                       for (i in seq_along(self$tracks)) {
                         if (is.null(self$tracks[[i]]$windows)) {
                           if(is.null(self$windows)) {stop("Global SeqPlot windows and at least one SeqTrack windows are NULL. All track windows must be set.")}
                           self$tracks[[i]]$windows <- self$windows
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

                       xscales <- unlist(lapply(self$tracks, function(t) {
                         lapply(seq_along(t$windows), function(i) {
                           c(start(t$windows)[i], end(t$windows)[i])
                         })
                       }), recursive = FALSE)

                       yscales <- list()

                       # axis expansion
                       #expandX <- if (!is.null(self$aesthetics$expandX)) self$aesthetics$expandX else c(0,0)
                       #expandY <- if (!is.null(self$aesthetics$expandY)) self$aesthetics$expandY else c(0,0)

                       #xscales <- lapply(xscales, expand_limits, expand = expandX)
                       #yscales <- lapply(yscales, expand_limits, expand = expandY)


                       for (track_idx in seq_along(self$tracks)) {
                         track <- self$tracks[[track_idx]]
                         disjoint <- isTRUE(track$aesthetics$disjointYScale)
                         n_windows <- length(track$windows)

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
                         grid_meta[[t]] <- vector("list", nWin)
                         x_left <- margins$left

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


                           track <- self$tracks[[t]]

                           grid_meta[[t]][[w]] <- list(
                             track  = t,
                             window = w,
                             full   = panel_full,
                             inner  = panel_full,
                             xscale = xscales[[panel_index]],
                             yscale = yscales[[panel_index]],
                             xScaleFactor = mcols(track$windows)$scale[[w]]
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
                     },

                     #' @description
                     #' Draw the grid backgrounds, track areas, window panels, and borders.
                     #' @return Renders the grid to the graphics device.
                     drawGrid = function() {
                       stopifnot(is.list(self$layout),
                                 !is.null(self$layout$panelBounds),
                                 !is.null(self$layout$trackBounds))

                       grid.newpage()
                       pushViewport(viewport(name = "root"))

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

                           # Draw window panel backgrounds inside the track
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

                       popViewport()
                     },

                     #' @description
                     #' Draw x- and y-axes for all tracks and windows, including ticks, labels,
                     #' and axis titles.
                     #' @return Renders axes to the graphics device.
                     drawAxes = function() {
                       panelBounds      <- self$layout$panelBounds
                       trackBounds <- self$layout$trackBounds
                       nTracks     <- length(panelBounds)

                       for (t in seq_len(nTracks)) {
                         trackAesthetics = modifyList(self$aesthetics, self$tracks[[t]]$aesthetics)

                         # Draw window axes over the panel borders
                         for (win in panelBounds[[t]]) {
                           xscale <- win$xscale
                           yscale <- win$yscale
                           p <- win$full

                           # x‐axis along bottom of this track
                           if (isTRUE(trackAesthetics$xAxisLine)) {
                             grid.lines(
                               x = unit(c(p$x0, p$x1), "npc"),
                               y = unit(c(p$y0, p$y0), "npc"),
                               gp = gpar(col = "#1C1B1A", lwd = 0.5)
                             )
                           }

                           # y‐axis along left of first panel
                           if (isTRUE(trackAesthetics$yAxisLine) &&
                               (isTRUE(trackAesthetics$yAxisPerWindow) || win$window == 1)) {
                             grid.lines(
                               x = unit(c(p$x0, p$x0), "npc"),
                               y = unit(c(p$y0, p$y1), "npc"),
                               gp = gpar(col = "#1C1B1A", lwd = 0.5)
                             )
                           }

                           # Draw x-axis ticks
                           if (isTRUE(trackAesthetics$xAxisTicks)) {
                             xbreaks <- xscale
                             xgrid   <- (xbreaks - xscale[1]) / diff(xscale)

                             sf <- if (!is.null(win$xScaleFactor)) win$xScaleFactor else 1e-6

                             unit_label <- if (abs(sf - 1e-6) < 1e-9) {
                               "Mb"
                             } else if (abs(sf - 1e-3) < 1e-9) {
                               "kb"
                             } else if (abs(sf - 1) < 1e-9) {
                               "bp"
                             } else {
                               paste0("×", signif(sf, 2))
                             }

                             xlabels <- paste0(format(round(xbreaks * sf), big.mark = ",", scientific = F), " ", unit_label)

                             for (i in seq_along(xgrid)) {
                               xpos <- win$full$x0 + xgrid[i] * (win$full$x1 - win$full$x0)
                               grid.lines(
                                 x = unit(c(xpos, xpos), "npc"),
                                 y = unit(c(win$full$y0, win$full$y0 - 0.005), "npc"),
                                 gp = gpar(col = "#1C1B1A", lwd = 0.5)
                               )

                               if (isTRUE(trackAesthetics$xAxisLabels)) {
                                 grid.text(
                                   #label = format(round(xbreaks[i] / 1000000), big.mark = ",", scientific = F),
                                   label = xlabels[i],
                                   x = unit(xpos, "npc"),
                                   y = unit(win$full$y0 - 0.015, "npc"),
                                   just = "top",
                                   rot = trackAesthetics$xAxisLabelRotation,
                                   hjust = trackAesthetics$xAxisLabelHorizontalJust,
                                   vjust = trackAesthetics$xAxisLabelVerticalJust,
                                   gp = gpar(cex = 0.6)
                                 )
                               }
                             }
                           }

                           # Draw y-axis ticks and labels
                           if (isTRUE(trackAesthetics$yAxisTicks) &&
                               (isTRUE(trackAesthetics$yAxisPerWindow) || win$window == 1)) {
                             ybreaks <- yscale
                             ygrid   <- (ybreaks - yscale[1]) / diff(yscale)
                             for (i in seq_along(ygrid)) {
                               ypos <- win$full$y0 + ygrid[i] * (win$full$y1 - win$full$y0)
                               grid.lines(
                                 x = unit(c(win$full$x0, win$full$x0 - 0.005), "npc"),
                                 y = unit(c(ypos, ypos), "npc"),
                                 gp = gpar(col = "#1C1B1A", lwd = 0.5)
                               )

                               if (isTRUE(trackAesthetics$yAxisLabels)) {
                                 grid.text(
                                   label = format(ybreaks[i], big.mark = ",", scientific = F),
                                   x = unit(win$full$x0 - 0.01, "npc"),
                                   y = unit(ypos, "npc"),
                                   just = "right",
                                   gp = gpar(cex = 0.6)
                                 )
                               }
                             }
                           }

                           # Draw x axis title
                           if (isTRUE(trackAesthetics$xAxisTitle)) {
                             x_label <- gsub("^chr", "", as.character(seqnames(self$tracks[[t]]$windows[win$window])))

                             grid.text(
                               label = x_label,
                               x = unit((p$x0 + p$x1) / 2, "npc"),
                               y = unit(p$y0 - 0.015, "npc"),
                               just = "top",
                               gp = gpar(cex = 0.6, fontface = "bold")
                             )
                           }

                           # Draw y axis title
                           if (isTRUE(trackAesthetics$yAxisTitle) && win$window == 1) {
                             y_label <- if (!is.null(trackAesthetics$yAxisTitleText)) {
                               trackAesthetics$yAxisTitleText
                             } else if (length(self$tracks[[t]]$elements) > 0 && !is.null(self$tracks[[t]]$elements[[1]]$yCol)) {
                               self$tracks[[t]]$elements[[1]]$yCol
                             } else {
                               ""
                             }

                             grid.text(
                               label = y_label,
                               x = unit(p$x0 - 0.03, "npc"),
                               y = unit((p$y0 + p$y1) / 2, "npc"),
                               just = "center",
                               rot = trackAesthetics$yAxisTitleRotation,
                               hjust = trackAesthetics$yAxisTitleHorizontalJust,
                               vjust = trackAesthetics$yAxisTitleVerticalJust,
                               gp = gpar(cex = 0.6, fontface = "bold")
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
                     #' Prepares and draws all components of the SeqPlot.
                     #' @return Renders elements to the graphics device.
                     plot = function() {
                       self$layoutGrid()
                       self$drawGrid()
                       self$drawAxes()
                       self$drawElements()
                       invisible(self)
                     }

                   ))


##' expand_limits
##' @description Expands axis limits around data
# expand_limits <- function(lims, expand) {
#   mult <- expand[1]
#   add  <- expand[2]
#   rng  <- diff(lims)
#   c(lims[1] - rng * mult - add,
#     lims[2] + rng * mult + add)
# }
#
