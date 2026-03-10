# Example 04: Position scales
# Exercises:
#   - seq_scale_genomic()  on SeqTrack (x-axis)
#   - seq_scale_continuous() on SeqTrack (y-axis, numeric data)
#   - seq_scale_discrete()  on SeqTrack (y-axis, categorical data)
#   - Reading scale_factor from mcols$scale
#
# Run with: devtools::load_all(); source("inst/examples/04_position_scales.R")

library(GenomicRanges)
devtools::load_all()


# ── seq_scale_genomic() ───────────────────────────────────────────────────────

# Default (1e-6 = Mb)
win_mb <- createGenomeWindows("chr1:1-100000000")
sg_mb  <- seq_scale_genomic(win_mb)
stopifnot(inherits(sg_mb, "SeqScaleGenomic"))
stopifnot(sg_mb$scale_factor == 1e-6)
message("04_position_scales: seq_scale_genomic default Mb OK")

# Explicit kb
sg_kb <- seq_scale_genomic(win_mb, scale_factor = 1e-3)
stopifnot(sg_kb$scale_factor == 1e-3)
message("04_position_scales: seq_scale_genomic explicit kb OK")

# Read from mcols$scale
win_custom         <- createGenomeWindows("chr1:1-100000000")
mcols(win_custom)$scale <- 1
sg_bp <- seq_scale_genomic(win_custom)
stopifnot(sg_bp$scale_factor == 1)
message("04_position_scales: seq_scale_genomic reads mcols$scale OK")

# Requires GRanges
tryCatch(
  seq_scale_genomic("not_a_granges"),
  error = function(e) message("04_position_scales: non-GRanges correctly errors")
)


# ── seq_scale_continuous() ────────────────────────────────────────────────────
sc_cont <- seq_scale_continuous(limits = c(-2, 2), n_breaks = 9)
stopifnot(inherits(sc_cont, "SeqScaleContinuous_Pos"))
stopifnot(sc_cont$limits == c(-2, 2))
stopifnot(sc_cont$n_breaks == 9)
message("04_position_scales: seq_scale_continuous OK")

# ── seq_scale_discrete() ─────────────────────────────────────────────────────
cell_types <- c("CD4_T", "CD8_T", "B_cell", "NK", "Mono")
sc_disc    <- seq_scale_discrete(
  levels = cell_types,
  labels = c("CD4+ T", "CD8+ T", "B", "NK", "Monocyte")
)
stopifnot(inherits(sc_disc, "SeqScaleDiscrete_Pos"))
stopifnot(identical(sc_disc$levels, cell_types))
stopifnot(length(sc_disc$labels) == 5L)
message("04_position_scales: seq_scale_discrete OK")

# ── Render: one track per scale type ─────────────────────────────────────────
win <- createGenomeWindows("chr1:1-100000000")

set.seed(7)
n   <- 50
snv <- GRanges("chr1", IRanges(sort(sample(1:1e8, n)), width = 1),
               score = rnorm(n, mean = 0, sd = 1))

tile_gr <- GRanges(
  "chr1",
  IRanges(seq(1, 100e6, by = 5e6), width = 4.5e6),
  cell_type = rep(cell_types, length.out = 20),
  color     = rep(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F"), length.out = 20)
)

plt <- SeqPlot(windows = win) %|%
  # Track 1: continuous y (z-score), genomic x (Mb, default)
  SeqTrack(scale_y = seq_scale_continuous(limits = c(-3, 3), n_breaks = 7),
           aesthetics = list(yAxisTitleText = "z-score")) %+%
    SeqPoint(snv, yCol = "score") %|%
  # Track 2: discrete y (cell types)
  SeqTrack(scale_y = seq_scale_discrete(levels = cell_types),
           aesthetics = list(yAxisTitleText = "Cell type")) %+%
    SeqTile(tile_gr, aes = aes(y = cell_type, fill = cell_type),
            scale = seq_scale_fill_discrete(
              values = setNames(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F"),
                                cell_types)))

grid::grid.newpage()
plt$layoutGrid()
plt$drawGrid()
plt$drawAxes()
plt$drawElements()
message("04_position_scales: rendered OK")
