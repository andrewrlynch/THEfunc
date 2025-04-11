#' Cytoband Table for hg38
#'
#' UCSC-style cytoband table for human genome (hg38).
#' Includes chromosomal location, name, and stain type.
#'
#' @format A data frame with 5 columns:
#' \describe{
#'   \item{chrom}{Chromosome name}
#'   \item{chromStart}{Start position}
#'   \item{chromEnd}{End position}
#'   \item{name}{Band name}
#'   \item{gieStain}{Staining pattern (e.g., gneg, gpos50)}
#' }
"cytoband_hg38"
