#' Load cytoband table
#'
#' This function loads a UCSC-style cytoband file.
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



#' Convert chromosome position to chromosome arm
#'
#' This function converts chromosome position to chromosome arm
#'
#' @param chr Chromosome
#' @param pos Position
#'
#' @return Character vector of chromosome arms
#' @export
positionToArm <- function(chr, pos, cyto = cytoband_hg38) {
  stopifnot(length(chr) == length(pos))  # enforce matching length

  # clean up chromosome labels
  cyto$chrom <- gsub("^chr", "", cyto$chrom, ignore.case = TRUE)
  chr <- gsub("chr|Chr|chromosome| ", "", chr)

  # for each query, find the matching cytoband
  res <- mapply(function(c, p) {
    hits <- cyto[cyto$chrom == c & cyto$chromStart <= p & cyto$chromEnd >= p, ]$name
    if (length(hits) == 0) return(NA_character_)
    substr(hits[1], 1, 1)  # first character: "p" or "q"
  }, chr, pos, USE.NAMES = FALSE)

  return(res)
}



#' Determine string-curve type from breakpoint strand
#'
#' Determines whether an SV link should be drawn as a C-shaped or S-shaped
#' cubic Bézier curve based on breakpoint strand orientation.
#'
#' Strand mapping (as requested):
#' \itemize{
#'   \item \code{+/+} → \code{"c"}
#'   \item \code{-/-} → \code{"c"}
#'   \item \code{+/-} → \code{"s"}
#'   \item \code{-/+} → \code{"s"}
#' }
#'
#' If either strand is \code{"*"} or not one of \code{"+"}/\code{"-"},
#' the function defaults to \code{default} (default: \code{"c"}).
#'
#' @param gr1 A \code{GRanges} for the first breakpoints.
#' @param gr2 A \code{GRanges} for the second breakpoints.
#' @param default Character scalar, \code{"c"} or \code{"s"}.
#'
#' @return A character vector of \code{"c"} or \code{"s"} with length \code{length(gr1)}.
#' @keywords internal
stringTypeFromStrand <- function(gr1, gr2, default = "c") {
  if (length(gr1) == 0) return(character(0))

  s1 <- as.character(GenomicRanges::strand(gr1))
  s2 <- as.character(GenomicRanges::strand(gr2))

  ok1 <- s1 %in% c("+", "-")
  ok2 <- s2 %in% c("+", "-")

  out <- rep(default, length(s1))
  out[ok1 & ok2 & (s1 != s2)] <- "s"
  out[ok1 & ok2 & (s1 == s2)] <- "c"
  out
}

#' Recycle/trim a vector to length \code{N}
#'
#' @param x Vector.
#' @param N Target length.
#' @param default Default value if \code{x} is \code{NULL} or length 0.
#'
#' @return Vector of length \code{N}.
#' @keywords internal
.recycled_to_N <- function(x, N, default) {
  if (N <= 0) return(x[0])
  if (is.null(x) || length(x) == 0) return(rep(default, N))
  if (length(x) == 1) return(rep(x, N))
  if (length(x) == N) return(x)
  rep(x, length.out = N)
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
  .Deprecated("createGenomeWindows",
    package = "THEfunc",
    msg = "createGenomeWindows is now the default function for creating genome windows. Please use createGenomeWindows instead of CreateSequenceWindows."
    )
  createGenomeWindows(regions = regions, padding = padding, genome = genome, add_chr = add_chr)
  }


#' createGenomeWindows
#' @name createGenomeWindows
#' @description
#' Generates a GRanges object of given sequence ranges
#' @param regions Vector of region strings such as "1:300-4000" or "chr2:20000-40000" or simply "chr3.
#' @param padding Number of base pairs around given positions to add. Automatically clipped to chromosome start positions. Default: 0
#' @param genome Genome of origin. Default: hg38
#' @param add_chr Whether to maintain "chr" convention on seqnames. Default: TRUE
#' @author Andrew Lynch
#' @export
createGenomeWindows <- function(regions, padding = 0, genome = "hg38", add_chr = TRUE) {
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
    gtf <- rtracklayer::import(gtf_path, format="gtf")
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
