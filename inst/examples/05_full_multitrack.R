# Example 05: Full multi-track plot — all new API features together
# Simulates a typical tumor sc-WGS overview:
#   Track 1 (Ideogram)  — chromosome bands
#   Track 2 (SeqRect)   — copy-number segments, colored by CN state
#   Track 3 (SeqPoint)  — SNV VAF, colored continuously by VAF (viridis)
#   Track 4 (SeqTile)   — single-cell CN heatmap, grouped + colored by clone
#   Track 5 (SeqBar)    — allele-specific coverage, colored by allele (discrete)
#
# Run with: devtools::load_all(); source("inst/examples/05_full_multitrack.R")

library(GenomicRanges)
devtools::load_all()

set.seed(99)
CHR  <- "chr8"
CLEN <- 145138636L   # hg38 chr8 length
win  <- CreateSequenceWindows(paste0(CHR, ":1-", CLEN))


# ── Simulate data ─────────────────────────────────────────────────────────────

## Copy-number segments (5 segments spanning chr8)
cn_starts <- c(1L,      20e6L, 55e6L, 90e6L, 120e6L)
cn_ends   <- c(19.9e6L, 54.9e6L, 89.9e6L, 119.9e6L, CLEN)
cn_vals   <- c(2L, 4L, 1L, 3L, 2L)
cn_colors <- c("2" = "#6BAED6", "4" = "#E6550D", "1" = "#74C476",
               "3" = "#9E9AC8", "0" = "#969696")
seg_gr <- GRanges(CHR, IRanges(cn_starts, cn_ends),
                  cn    = cn_vals,
                  color = cn_colors[as.character(cn_vals)])

## SNVs (60 random positions, VAF from Beta)
n_snv <- 60
snv_gr <- GRanges(
  CHR,
  IRanges(sort(sample(1:CLEN, n_snv)), width = 1),
  vaf = round(rbeta(n_snv, 2, 4), 3)
)

## Single-cell CN heatmap (4 clones × 30 bins)
clones     <- c("Clone_A", "Clone_B", "Clone_C", "Clone_D")
clone_cols <- c(Clone_A = "#4E79A7", Clone_B = "#F28E2B",
                Clone_C = "#E15759", Clone_D = "#59A14F")
n_bins     <- 30
bin_starts <- seq(1L, CLEN - 5e6L, length.out = n_bins)
tile_gr <- GRanges(
  CHR,
  IRanges(rep(as.integer(bin_starts), times = length(clones)),
          width = as.integer(CLEN / n_bins * 0.9)),
  clone = rep(clones, each = n_bins)
)

## Allele-specific bars (major / minor allele depth per segment)
bar_gr <- GRanges(
  CHR,
  IRanges(cn_starts, cn_ends),
  major = as.numeric(cn_vals) * 0.6 + rnorm(5, 0, 0.05),
  allele = rep(c("major", "minor"), length.out = 5),
  color  = rep(c("#4E79A7", "#F28E2B"), length.out = 5)
)


# ── Build plot ────────────────────────────────────────────────────────────────

plt <- SeqPlot(windows = win) %|%

  # Track 1: Ideogram (uses package cytoband data)
  SeqTrack() %+%
    SeqIdeogram(cytoband_hg38) %|%

  # Track 2: Copy-number segments, colored by CN (pre-assigned in color mcol)
  SeqTrack(scale_y = seq_scale_continuous(limits = c(0, 5), n_breaks = 6),
           aesthetics = list(yAxisTitleText = "CN")) %+%
    SeqRect(seg_gr, yCol = "cn") %|%

  # Track 3: SNV VAF, continuously colored by VAF via aes()
  SeqTrack(scale_y = seq_scale_continuous(limits = c(0, 1), n_breaks = 5),
           aesthetics = list(yAxisTitleText = "VAF")) %+%
    SeqPoint(
      snv_gr,
      yCol  = "vaf",
      aes   = aes(color = vaf),
      scale = seq_scale_color_continuous(palette = "viridis", limits = c(0, 1))
    ) %|%

  # Track 4: Single-cell clone heatmap — discrete fill via aes()
  SeqTrack(scale_y = seq_scale_discrete(levels = clones),
           aesthetics = list(yAxisTitleText = "Clone")) %+%
    SeqTile(
      tile_gr,
      aes   = aes(y = clone, fill = clone),
      scale = seq_scale_fill_discrete(values = clone_cols)
    ) %|%

  # Track 5: Allele bar (color pre-set in gr mcol)
  SeqTrack(scale_y = seq_scale_continuous(limits = c(0, 3), n_breaks = 4),
           aesthetics = list(yAxisTitleText = "Depth")) %+%
    SeqBar(bar_gr, yCol = "major")


# ── Structural assertions ─────────────────────────────────────────────────────
stopifnot(length(plt$tracks) == 5L)
stopifnot(all(vapply(plt$tracks, function(t) length(t$elements) == 1L, logical(1))))
message("05_full_multitrack: structure checks passed")

# Clone tile groups correct
tile_el <- plt$tracks[[4]]$elements[[1]]
stopifnot(identical(tile_el$groups, clones))
message("05_full_multitrack: clone grouping OK")

# VAF point colors are hex strings
pt_el <- plt$tracks[[3]]$elements[[1]]
stopifnot(all(grepl("^#[0-9A-Fa-f]{6}$", pt_el$aesthetics$color)))
message("05_full_multitrack: continuous color OK")


# ── Render ────────────────────────────────────────────────────────────────────
grid::grid.newpage()
plt$layoutGrid(trackHeights = c(0.3, 1, 1.5, 1.5, 1))
plt$drawGrid()
plt$drawAxes()
plt$drawElements()
message("05_full_multitrack: rendered OK")
