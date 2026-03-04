# Example 02: SeqTile with aes(y = ..., fill = ...) — discrete fill
# Exercises:
#   - aes(y = cell_type)  → groups / y-index from a character column
#   - aes(fill = cell_type) + seq_scale_fill_discrete(values = ...) → named colors
#   - seq_scale_discrete() on the track y-axis
#
# Run with: devtools::load_all(); source("inst/examples/02_aes_discrete_tile.R")

library(GenomicRanges)
devtools::load_all()

win <- CreateSequenceWindows("chr1:1-100000000")

set.seed(1)
n_cells  <- 6
n_bins   <- 40
cell_ids <- paste0("Cell_", LETTERS[1:n_cells])

# Simulate copy-number tiles: each cell × genomic bin
tile_gr <- GRanges(
  seqnames  = "chr1",
  ranges    = IRanges(
    start = rep(seq(1, 95e6, length.out = n_bins), times = n_cells),
    width = 2.5e6
  ),
  cell_type = rep(cell_ids, each = n_bins)
)

# ── Build tile element ────────────────────────────────────────────────────────
cell_colors <- c(
  Cell_A = "#4E79A7", Cell_B = "#F28E2B", Cell_C = "#E15759",
  Cell_D = "#76B7B2", Cell_E = "#59A14F", Cell_F = "#EDC948"
)

tile <- SeqTile(
  tile_gr,
  aes   = aes(y = cell_type, fill = cell_type),
  scale = seq_scale_fill_discrete(values = cell_colors)
)

# Structural checks
stopifnot(identical(tile$groups, cell_ids))
stopifnot(all(tile$y == rep(seq_len(n_cells), each = n_bins)))
stopifnot(tile$groupCol == "cell_type")
message("02_aes_discrete_tile: grouping checks passed")

# Color checks — Cell_A tiles should all be #4E79A7
cell_a_idx <- which(as.character(mcols(tile$gr)$cell_type) == "Cell_A")
stopifnot(all(mcols(tile$gr)$color[cell_a_idx] == "#4E79A7"))
message("02_aes_discrete_tile: color checks passed")

# ── Inferred y scale ──────────────────────────────────────────────────────────
inferred <- tile$.infer_scale_y()
stopifnot(inherits(inferred, "SeqScaleDiscrete_Pos"))
stopifnot(identical(inferred$levels, cell_ids))
message("02_aes_discrete_tile: inferred y scale OK")

# ── Render ────────────────────────────────────────────────────────────────────
plt <- SeqPlot(windows = win) %|%
  SeqTrack(scale_y = inferred) %+% tile

grid::grid.newpage()
plt$layoutGrid()
plt$drawGrid()
plt$drawAxes()
plt$drawElements()
message("02_aes_discrete_tile: rendered OK")
