# Example 01: Operator pipeline (%|% and %+%)
# Verifies that SeqPlot() %|% SeqTrack() %+% SeqElement() builds the correct
# object structure, and that the old $new() style now correctly errors.
#
# Run with: devtools::load_all(); source("inst/examples/01_operator_pipeline.R")

library(GenomicRanges)
devtools::load_all()

win <- createGenomeWindows("chr7:1-159345973")

# Synthetic SNV positions with a numeric score column
snv <- GRanges(
  seqnames = "chr7",
  ranges   = IRanges(
    start = c(5e6, 15e6, 40e6, 72e6, 110e6, 140e6),
    width = 1
  ),
  vaf = c(0.12, 0.48, 0.33, 0.51, 0.27, 0.44)
)

# Copy-number segments
segs <- GRanges(
  seqnames = "chr7",
  ranges   = IRanges(
    start = c(1,    30e6, 65e6, 100e6),
    end   = c(29e6, 64e6, 99e6, 159e6)
  ),
  cn    = c(2L, 4L, 1L, 3L),
  color = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2")
)


# ── Build the plot with the new operator API ──────────────────────────────────
plt <- SeqPlot(windows = win) %|%
  SeqTrack(scale_y = seq_scale_continuous(limits = c(0, 6), n_breaks = 4)) %+%
    SeqRect(segs, yCol = "cn") %|%
  SeqTrack(scale_y = seq_scale_continuous(limits = c(0, 1), n_breaks = 3)) %+%
    SeqPoint(snv, yCol = "vaf")

stopifnot(length(plt$tracks) == 2L)
stopifnot(length(plt$tracks[[1]]$elements) == 1L)
stopifnot(length(plt$tracks[[2]]$elements) == 1L)
stopifnot(inherits(plt$tracks[[1]]$elements[[1]], "SeqRect"))
stopifnot(inherits(plt$tracks[[2]]$elements[[1]], "SeqPoint"))
message("01_operator_pipeline: structure checks passed")

# ── Verify that $new() on a wrapped class now errors ─────────────────────────
tryCatch(
  SeqPoint$new(snv, yCol = "vaf"),
  error = function(e) message("01_operator_pipeline: $new() correctly errors — ", conditionMessage(e))
)

# ── Render ────────────────────────────────────────────────────────────────────
grid::grid.newpage()
plt$layoutGrid()
plt$drawGrid()
plt$drawAxes()
plt$drawElements()
message("01_operator_pipeline: rendered OK")
