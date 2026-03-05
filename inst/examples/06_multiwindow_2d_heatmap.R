# Two-dimensional genomic heatmaps with SeqTile
#
# This example demonstrates the new x/y axes feature enabling:
#   - 2D genomic heatmaps (genomic x and y axes)
#   - Multi-window layouts (showing multiple disjoint regions)
#   - Contact-style plots for Hi-C, SV co-occurrence, etc.
#
# Key concepts:
#   - SeqTile$new(x=gr_x, y=gr_y) creates a 2D tile element
#   - SeqTrack(windows=x_windows, y_windows=y_windows) sets up genomic axes
#   - All constructors now use x=/x1=/x2= (gr=/gr1=/gr2= still work for backwards compat)

library(GenomicRanges)
library(IRanges)
library(grid)

# ══════════════════════════════════════════════════════════════════════════════
# Example 1: Simple 2D genomic heatmap (chr1 vs chr2)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Example 1: 2D Genomic Heatmap (chr1 × chr2 contact map) ===\n")

set.seed(42)

# Simulate chr1 × chr2 contact map (100 tiles)
n_contacts <- 100

# X-axis: contact points on chr1 (1-50 Mb)
x_starts <- sort(sample(seq(1e6, 50e6, by=500e3), n_contacts, replace=TRUE))
x_gr <- GRanges(
  seqnames = "chr1",
  ranges   = IRanges(start = x_starts, width = 1e6)
)

# Y-axis: contact points on chr2 (5-40 Mb)
y_starts <- sort(sample(seq(5e6, 40e6, by=500e3), n_contacts, replace=TRUE))
y_gr <- GRanges(
  seqnames = "chr2",
  ranges   = IRanges(start = y_starts, width = 1e6)
)

# Contact strength colorized (blue = weak, red = strong)
strength <- runif(n_contacts, 0, 1)
cols <- colorRampPalette(c("#E8F4FF", "#4385BE", "#205EA6"))(100)[pmax(1, round(strength * 99) + 1)]
mcols(x_gr)$color <- cols

# Create 2D tile element: both axes are genomic
tile_2d <- SeqTile(
  x = x_gr,           # X-axis: genomic positions on chr1
  y = y_gr,           # Y-axis: genomic positions on chr2
  aesthetics = list(border = NA, lwd = 0.1)
)

# Define viewing windows
x_window <- GRanges("chr1", IRanges(1e6, 51e6))
mcols(x_window)$scale <- 1e-6  # Display in Mb

y_window <- GRanges("chr2", IRanges(4e6, 41e6))
mcols(y_window)$scale <- 1e-6

# Create track with both x-windows (genomic) and y-windows (genomic)
track_2d <- SeqTrack(
  windows   = x_window,
  y_windows = y_window,
  aesthetics = list(
    xAxisTitle = TRUE,
    xAxisTitleText = "chr1 (Mb)",
    yAxisTitle = TRUE,
    yAxisTitleText = "chr2 (Mb)"
  )
)

# Add the 2D tile to the track
track_2d$addElement(tile_2d)

# Create plot
sp1 <- SeqPlot(
  tracks = list(track_2d),
  windows = x_window,
  aesthetics = list(
    margins = list(top = 0.10, right = 0.10, bottom = 0.12, left = 0.12),
    trackHeights = 1,
    xAxisLabels = TRUE,
    yAxisLabels = TRUE
  )
)

# Render to screen
sp1$plot()
grid.text(
  "Example 1: 2D Genomic Heatmap (chr1 × chr2)",
  x = 0.5, y = 0.98, just = "top",
  gp = gpar(fontsize = 14, fontface = "bold")
)

# ══════════════════════════════════════════════════════════════════════════════
# Example 2: Multi-window layout (two disjoint x-regions, same y-region)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Example 2: Multi-window (chr3 regions A & B vs chr4) ===\n")

set.seed(99)

# Two separate regions on chr3: 1-15 Mb and 30-45 Mb
n_region_a <- 50
n_region_b <- 50

# Region A
x_a_starts <- sort(sample(seq(1e6, 15e6, by=500e3), n_region_a, replace=TRUE))
x_a <- GRanges("chr3", IRanges(x_a_starts, width=1e6))

# Region B
x_b_starts <- sort(sample(seq(30e6, 45e6, by=500e3), n_region_b, replace=TRUE))
x_b <- GRanges("chr3", IRanges(x_b_starts, width=1e6))

x_all <- c(x_a, x_b)

# Y-axis: chr4 (2-20 Mb)
y_starts <- sort(sample(seq(2e6, 20e6, by=500e3), n_region_a + n_region_b, replace=TRUE))
y_all <- GRanges("chr4", IRanges(y_starts, width=1e6))

# Color map
vals <- runif(n_region_a + n_region_b, 0, 1)
cols_multi <- colorRampPalette(c("white", "#D14D41", "#A02C6D"))(100)[pmax(1, round(vals * 99) + 1)]
mcols(x_all)$color <- cols_multi

tile_multi <- SeqTile(
  x = x_all,
  y = y_all,
  aesthetics = list(border = NA)
)

# Two x-windows (Region A and Region B)
x_windows <- GRanges(
  seqnames = c("chr3", "chr3"),
  ranges = IRanges(c(1e6, 30e6), end = c(15e6, 45e6))
)
mcols(x_windows)$scale <- rep(1e-6, 2)

# Single y-window
y_window_multi <- GRanges("chr4", IRanges(1e6, 21e6))
mcols(y_window_multi)$scale <- 1e-6

track_multi <- SeqTrack(
  windows = x_windows,
  y_windows = y_window_multi,
  aesthetics = list(
    xAxisTitle = TRUE,
    xAxisTitleText = "chr3 (Mb)",
    yAxisTitle = TRUE,
    yAxisTitleText = "chr4 (Mb)"
  )
)

track_multi$addElement(tile_multi)

sp2 <- SeqPlot(
  tracks = list(track_multi),
  windows = x_windows,
  aesthetics = list(
    margins = list(top = 0.10, right = 0.10, bottom = 0.12, left = 0.12),
    windowGaps = 0.03,
    trackHeights = 1,
    xAxisLabels = TRUE,
    yAxisLabels = TRUE
  )
)

sp2$plot()
grid.text(
  "Example 2: Multi-window Heatmap (chr3 regions A & B vs chr4)",
  x = 0.5, y = 0.98, just = "top",
  gp = gpar(fontsize = 14, fontface = "bold")
)

# ══════════════════════════════════════════════════════════════════════════════
# Example 3: Multi-window on both axes (2×2 grid view)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Example 3: Multi-window on both axes (2×2 grid) ===\n")

set.seed(77)

# Two x-windows on chr5: 1-12 Mb and 20-32 Mb
n_x_a <- 40
n_x_b <- 40
n_total <- n_x_a + n_x_b

x_grid_a <- sort(sample(seq(1e6, 12e6, by=500e3), n_x_a, replace=TRUE))
x_grid_a_gr <- GRanges("chr5", IRanges(x_grid_a, width=1e6))

x_grid_b <- sort(sample(seq(20e6, 32e6, by=500e3), n_x_b, replace=TRUE))
x_grid_b_gr <- GRanges("chr5", IRanges(x_grid_b, width=1e6))

x_grid_all <- c(x_grid_a_gr, x_grid_b_gr)

# Two y-windows on chr6: 1-10 Mb and 15-25 Mb
# Create exactly n_total y-axis contact points total
y_grid_a <- sort(sample(seq(1e6, 10e6, by=500e3), n_x_a, replace=TRUE))
y_grid_b <- sort(sample(seq(15e6, 25e6, by=500e3), n_x_b, replace=TRUE))
y_grid_all_vec <- c(y_grid_a, y_grid_b)
y_grid_all <- GRanges("chr6", IRanges(y_grid_all_vec, width=1e6))

# Colorize with proper length
vals_grid <- runif(n_total, 0, 1)
cols_grid <- colorRampPalette(c("#EEF3FF", "#66A0C8", "#205EA6"))(100)[pmax(1, round(vals_grid * 99) + 1)]
mcols(x_grid_all)$color <- cols_grid

tile_grid <- SeqTile(
  x = x_grid_all,
  y = y_grid_all,
  aesthetics = list(border = NA)
)

# Two x-windows
x_wins_grid <- GRanges(
  seqnames = c("chr5", "chr5"),
  ranges = IRanges(c(1e6, 20e6), end = c(12e6, 32e6))
)

# Two y-windows (note: order corresponds to how y_grid_mixed is structured)
y_wins_grid <- GRanges(
  seqnames = c("chr6", "chr7", "chr8"),
  ranges = IRanges(c(1e6, 15e6, 25e6), end = c(10e6, 25e6, 35e6))
)

track_grid <- SeqTrack(
  windows = x_wins_grid,
  y_windows = y_wins_grid,
  aesthetics = list(
    xAxisTitle = TRUE,
    yAxisTitle = TRUE,
    yWindowGap = 0.02
  )
)

track_grid$addElement(tile_grid)

sp3 <- SeqPlot(
  tracks = list(track_grid),
  windows = x_wins_grid,
  aesthetics = list(
    margins = list(top = 0.10, right = 0.10, bottom = 0.12, left = 0.12),
    windowGaps = 0.02,
    trackHeights = 1,
    xAxisLabels = TRUE,
    yAxisLabels = TRUE,
    windowBorder = TRUE,
    trackBorder = FALSE,
    yAxisCap = "exact"
  )
)

sp3$plot()
grid.text(
  "Example 3: Multi-window 2×2 Grid (chr5 × chr6)",
  x = 0.5, y = 0.98, just = "top",
  gp = gpar(fontsize = 14, fontface = "bold")
)

# ══════════════════════════════════════════════════════════════════════════════
# Example 4: Backwards compatibility - 1D SeqTile still works
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Example 4: 1D SeqTile (categorical y-axis, backwards compatible) ===\n")

set.seed(55)

n_categories <- 60
x_cat_starts <- sort(sample(seq(1e6, 30e6, by=500e3), n_categories, replace=TRUE))
sample_types <- sample(c("TypeA", "TypeB", "TypeC"), n_categories, replace=TRUE)

tile_colors_cat <- c(
  TypeA = "#4385BE",
  TypeB = "#879A39",
  TypeC = "#D14D41"
)[sample_types]

gr_cat <- GRanges(
  "chr7",
  IRanges(x_cat_starts, width=1e6),
  type = sample_types,
  color = tile_colors_cat
)

# 1D mode: only x is genomic, y is categorical
tile_1d <- SeqTile(
  x = gr_cat,           # Genomic x-axis
  groupCol = "type",    # Categorical y-axis
  aesthetics = list(border = NA)
)

track_1d <- SeqTrack(
  windows = GRanges("chr7", IRanges(1e6, 31e6)),
  aesthetics = list(
    xAxisTitle = TRUE,
    xAxisTitleText = "chr7 (Mb)",
    yAxisTitle = TRUE,
    yAxisTitleText = "Sample Type"
  )
)

track_1d$addElement(tile_1d)

sp4 <- SeqPlot(
  tracks = list(track_1d),
  windows = GRanges("chr7", IRanges(1e6, 31e6)),
  aesthetics = list(
    margins = list(top = 0.10, right = 0.10, bottom = 0.12, left = 0.12),
    trackHeights = 1,
    xAxisLabels = TRUE,
    yAxisLabels = TRUE
  )
)

sp4$plot()
grid.text(
  "Example 4: 1D SeqTile (categorical y-axis) - Backwards Compat",
  x = 0.5, y = 0.98, just = "top",
  gp = gpar(fontsize = 14, fontface = "bold")
)

cat("\n✓ All examples completed!\n")
