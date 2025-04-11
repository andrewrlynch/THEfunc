# -----------------------------------------------------------------------------
# Generate random grange data
# -----------------------------------------------------------------------------
make_random_hg38_granges <- function(n_intervals = 1000, min_width = 50, max_width = 1000) {
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
