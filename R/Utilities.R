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
