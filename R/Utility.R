#' makeRandomHg38Grange
#' @description Generate a random hg38 grange object
#' @export
makeRandomHg38Grange <- function(n_intervals = 1000, min_width = 50, max_width = 1000) {
  # Standard hg38 chromosome lengths (chr1-22, X, Y)
  chr_lengths <- c(
    "chr1"  = 248956422,
    "chr2"  = 242193529,
    "chr3"  = 198295559,
    "chr4"  = 190214555,
    "chr5"  = 181538259,
    "chr6"  = 170805979,
    "chr7"  = 159345973,
    "chr8"  = 145138636,
    "chr9"  = 138394717,
    "chr10" = 133797422,
    "chr11" = 135086622,
    "chr12" = 133275309,
    "chr13" = 114364328,
    "chr14" = 107043718,
    "chr15" = 101991189,
    "chr16" = 90338345,
    "chr17" = 83257441,
    "chr18" = 80373285,
    "chr19" = 58617616,
    "chr20" = 64444167,
    "chr21" = 46709983,
    "chr22" = 50818468,
    "chrX"  = 156040895,
    "chrY"  = 57227415
  )


  # Sample chromosomes with probability proportional to their lengths
  chrs <- sample(names(chr_lengths), size = n_intervals, replace = TRUE, prob = chr_lengths)

  # Preallocate vectors for starts, widths, and ends
  starts <- numeric(n_intervals)
  widths <- numeric(n_intervals)
  ends <- numeric(n_intervals)

  for (i in seq_len(n_intervals)) {
    chr <- chrs[i]
    chr_len <- chr_lengths[chr]

    # Randomly choose a width within the specified range
    width_candidate <- sample(min_width:max_width, 1)
    # Ensure that the chosen width does not exceed the chromosome length
    if (width_candidate > chr_len) {
      width_candidate <- chr_len
    }

    # Maximum allowed start such that the interval does not exceed chromosome length
    max_start <- chr_len - width_candidate + 1
    # Randomly choose a start position in the allowed range
    start_val <- sample(seq_len(max_start), 1)
    end_val <- start_val + width_candidate - 1

    starts[i] <- start_val
    widths[i] <- width_candidate
    ends[i] <- end_val
  }

  # Create GRanges object with the random intervals
  gr <- GRanges(seqnames = chrs,
                ranges = IRanges(start = starts, end = ends))

  # Add a random score (between 0 and 10) to the metadata columns
  mcols(gr)$score <- runif(n_intervals, 0, 10)

  return(gr)
}



#' AnnotateSvWithSegmentCn
#' @description Annotate an SV vcf with total copy number at breakend positions.
#' @export
AnnotateSvWithSegmentCn <- function(sv_gr, cn_gr, slop = 5000, cn_col = "copy_number") {

  sv_ranges <- StructuralVariantAnnotation::breakpointRanges(sv_gr)
  sv_pairs <- StructuralVariantAnnotation::breakpointgr2pairs(sv_ranges)

  # Helper function for annotating one end
  annotate_one_end <- function(bkpt_gr, cn_gr, slop, cn_col) {
    sv_slop <- suppressWarnings(resize(bkpt_gr, width = 1, fix = "start"))
    sv_slop <- flank(sv_slop, width = slop, both = TRUE)

    hits <- findOverlaps(sv_slop, cn_gr)
    cn_vals <- rep(NA_real_, length(bkpt_gr))

    for (i in seq_along(bkpt_gr)) {
      print(i)
      olap_idx <- subjectHits(hits)[queryHits(hits) == i]
      if (length(olap_idx) == 0) next

      overlapping_cns <- cn_gr[olap_idx]
      cn_values <- mcols(overlapping_cns)[[cn_col]]

      strand_val <- as.character(strand(bkpt_gr[i]))

      if (length(overlapping_cns) == 1) {
        cn_vals[i] <- cn_values
      } else if (strand_val == "+") {
        cn_vals[i] <- cn_values[which.min(start(overlapping_cns))]
      } else if (strand_val == "-") {
        cn_vals[i] <- cn_values[which.max(end(overlapping_cns))]
      } else {
        cn_vals[i] <- cn_values[1]  # fallback for strand == "*"
      }
    }

    return(cn_vals)
  }

  # Annotate first and second breakends
  cn_first <- annotate_one_end(sv_pairs@first, cn_gr, slop, cn_col)
  cn_second <- annotate_one_end(sv_pairs@second, cn_gr, slop, cn_col)

  # Add annotations to metadata
  mcols(sv_pairs)$cn_first <- cn_first
  mcols(sv_pairs)$cn_second <- cn_second

  return(sv_pairs)
}

#' CreateSequenceWindows
#' @description Create a GRanges object with specified ranges using e.g., chr1:100-200
#' @export
CreateSequenceWindows <- function(regions, padding = 0, genome = "hg38", add_chr = FALSE) {
  # Ensure input is character
  regions <- as.character(regions)

  # Pattern to extract chrom, start, and end
  parsed <- strcapture(
    pattern = "^([^:]+):([0-9,]+)-([0-9,]+)$",
    x = regions,
    proto = list(chrom = character(), start = character(), end = character())
  )

  # Sanity check
  if (nrow(parsed) != length(regions)) {
    stop("Could not parse one or more regions. Make sure they are like 'chr1:100-200' or '1:100-200'.")
  }

  # Remove commas and convert to integers
  parsed$start <- as.integer(gsub(",", "", parsed$start)) - padding
  parsed$end   <- as.integer(gsub(",", "", parsed$end)) + padding
  parsed$start[parsed$start < 1] <- 1

  chrom <- if (add_chr) paste0("chr", parsed$chrom) else parsed$chrom

  GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = parsed$start, end = parsed$end)
  ) |> sort()
}


