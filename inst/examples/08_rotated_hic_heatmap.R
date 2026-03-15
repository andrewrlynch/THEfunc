# Example 08: Rotated Hi-C Heatmap Visualization
#
# This example demonstrates practical visualization of symmetric genomic interaction
# matrices (Hi-C contact maps) using SeqTile's rotated styles.
#
# The rotated styles enable space-efficient visualization of symmetric data:
#   - style="full": Traditional rectangular heatmap (shows full matrix)
#   - style="diagonal": Lower diagonal only (removes redundant upper triangle)
#   - style="triangle": 45° rotated forming triangle (compact representation)
#   - style="rectangle": 45° rotated with window expansion (shows off-diagonal contacts)
#
# Run with: devtools::load_all(); source("inst/examples/08_rotated_hic_heatmap.R")

library(GenomicRanges)
library(IRanges)
library(grid)

devtools::load_all()

cat("\n=== Hi-C Heatmap Visualization with Rotated Styles ===\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 1. Create Realistic Hi-C Contact Data
# ──────────────────────────────────────────────────────────────────────────────

cat("1. Creating Hi-C contact matrix...\n")

set.seed(123)

# Genomic region: chr5:40Mb-50Mb (10 Mb region)
# Bin size: 100 kb
region_start <- 40e6
region_end <- 50e6
bin_size <- 100000
n_bins <- (region_end - region_start) / bin_size

# Create bin boundaries
bin_starts <- seq(region_start, region_end - bin_size, by = bin_size)
bin_ends <- bin_starts + bin_size

cat("  Region: chr5:", region_start/1e6, "-", region_end/1e6, " Mb\n")
cat("  Bin size:", bin_size/1e3, " kb\n")
cat("  Number of bins:", n_bins, "\n\n")

# Generate contact frequencies with realistic patterns:
# - Diagonal contacts are strong (same bin is most frequent)
# - Off-diagonal strength decays with distance
# - Add some noise
generate_hic_matrix <- function(n_bins, decay_rate = 0.2) {
  contacts <- list()

  for (i in 1:n_bins) {
    for (j in i:n_bins) {
      # Distance decay: exponential falloff
      distance <- j - i
      base_strength <- exp(-distance * decay_rate)

      # Add stochasticity
      strength <- base_strength * rlnorm(1, meanlog = 0, sdlog = 0.3)
      strength <- pmax(0.01, strength)  # Minimum contact strength

      # Store contact
      contacts[[length(contacts) + 1]] <- list(
        bin_i = i,
        bin_j = j,
        strength = strength
      )

      # Add symmetric contact (if not diagonal)
      if (i != j) {
        contacts[[length(contacts) + 1]] <- list(
          bin_i = j,
          bin_j = i,
          strength = strength
        )
      }
    }
  }

  data.frame(
    bin_i = sapply(contacts, function(x) x$bin_i),
    bin_j = sapply(contacts, function(x) x$bin_j),
    strength = sapply(contacts, function(x) x$strength)
  )
}

hic_matrix <- generate_hic_matrix(n_bins, decay_rate = 0.25)

cat("Generated", nrow(hic_matrix), "contact pairs\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Create Color-Mapped Tiles
# ──────────────────────────────────────────────────────────────────────────────

cat("2. Mapping contacts to colors...\n")

# Create color palette (blue=weak, red=strong)
color_palette <- colorRampPalette(RColorBrewer::brewer.pal(5,"RdPu"))(256)

# Map strength values to colors
strength_scaled <- (hic_matrix$strength - min(hic_matrix$strength)) /
  (max(hic_matrix$strength) - min(hic_matrix$strength))
hic_matrix$color <- color_palette[pmax(1, round(strength_scaled * 255) + 1)]

# Create GRanges for x and y axes
x_gr <- GRanges(
  seqnames = "chr5",
  ranges = IRanges(
    start = bin_starts[hic_matrix$bin_i],
    width = bin_size
  ),
  strength = hic_matrix$strength,
  color = hic_matrix$color
)

y_gr <- GRanges(
  seqnames = "chr5",
  ranges = IRanges(
    start = bin_starts[hic_matrix$bin_j],
    width = bin_size
  ),
  strength = hic_matrix$strength,
  color = hic_matrix$color
)

cat("  Color palette: blue (weak) → red (strong)\n")
cat("  X-axis: genomic positions\n")
cat("  Y-axis: genomic positions\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 3. Create Visualization Window
# ──────────────────────────────────────────────────────────────────────────────

window <- GRanges("chr5", IRanges(region_start, region_end))
mcols(window)$scale <- 1e-6  # Display in Mb

cat("3. Creating visualization window...\n")
cat("  Window: chr5:", start(window)/1e6, "-", end(window)/1e6, " Mb\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Create and Render All 4 Styles
# ──────────────────────────────────────────────────────────────────────────────

cat("4. Creating plots for all 4 styles...\n\n")

# Style 1: Full heatmap
cat("   Style 1: FULL (Traditional rectangular heatmap)\n")
cat("   - Shows complete symmetric matrix\n")
cat("   - Useful for detailed analysis\n")
cat("   - Requires more space (redundant upper triangle)\n\n")

tile_full <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "full",
  aesthetics = list(border = NA, lwd = 0.05)
)

plt_full <- SeqPlot(windows = window) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
  tile_full

cat("   Rendering full style plot...\n")
plt_full$layoutGrid()
plt_full$drawGrid()
plt_full$drawAxes()
plt_full$drawElements()

# Style 2: Diagonal only
cat("\n   Style 2: DIAGONAL (Lower diagonal only)\n")
cat("   - Removes redundant upper triangle\n")
cat("   - Still uses rectangular tiles\n")
cat("   - More compact than full\n\n")

tile_diag <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "diagonal",
  aesthetics = list(border = NA, lwd = 0.05)
)

plt_diag <- SeqPlot(windows = window) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
  tile_diag

cat("   Rendering diagonal style plot...\n")
plt_diag$layoutGrid()
plt_diag$drawGrid()
plt_diag$drawAxes()
plt_diag$drawElements()

# Style 3: Triangle (45° rotated)
cat("\n   Style 3: TRIANGLE (45° rotated, compact)\n")
cat("   - Lower diagonal rotated 45°\n")
cat("   - Forms triangular shape\n")
cat("   - Most space-efficient\n")
cat("   - Standard for Hi-C visualization\n\n")

tile_triangle <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "triangle",
  yCoordType = "distance",
  yDistMax = 1000000,
  aesthetics = list(border = NA, lwd = 0.1)
)

plt_triangle <- SeqPlot(windows = window, aesthetics = list(margins = list(top = 0.05, bottom = 0.1, left = 0.1, right = 0.05))) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA), windows = createGenomeWindows("chr5:44000000-48000000")) %+%
  tile_triangle

cat("   Rendering triangle style plot...\n")
plt_triangle$plot()

# Style 4: Rectangle (45° rotated with expansion)
cat("\n   Style 4: RECTANGLE (45° rotated with window expansion)\n")
cat("   - Rotated lower diagonal in rectangular frame\n")
cat("   - Window expansion shows off-diagonal contacts\n")
cat("   - Useful for studying TAD (topologically associating domain) boundaries\n\n")

maxDist <- 2e6  # 2 Mb expansion on each side
tile_rect <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "rectangle",
  yCoordType = "distance",
  yDistMax = 1000000,
  aesthetics = list(border = NA, lwd = 0.05)
)

plt_rect <- SeqPlot(aesthetics = list(margins = list(top = 0.1, bottom = 0.2, left = 0.1, right = 0.1))) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA), windows = createGenomeWindows(c("chr5:42000000-45000000"))) %+%
  tile_rect

cat("   Window expansion: ±", maxDist/1e6, " Mb\n")
cat("   Rendering rectangle style plot...\n")
plt_rect$layoutGrid()
plt_rect$drawGrid()
plt_rect$drawAxes()
plt_rect$drawElements()

# ──────────────────────────────────────────────────────────────────────────────
# 5. Summary and Comparison
# ──────────────────────────────────────────────────────────────────────────────

cat("\n=== Summary: Choosing the Right Style ===\n\n")

cat("FULL:\n")
cat("  ✓ Shows complete picture\n")
cat("  ✓ Good for presentations/publications\n")
cat("  ✗ Wastes space on redundant data\n\n")

cat("DIAGONAL:\n")
cat("  ✓ Removes redundancy\n")
cat("  ✓ Traditional rectangular tiles\n")
cat("  ✗ Still uses considerable space\n\n")

cat("TRIANGLE:\n")
cat("  ✓ Most compact (45° rotation)\n")
cat("  ✓ Standard in bioinformatics\n")
cat("  ✓ Shows local TAD structure well\n")
cat("  ✗ Less intuitive for unfamiliar viewers\n\n")

cat("RECTANGLE:\n")
cat("  ✓ Shows off-diagonal region\n")
cat("  ✓ Window expansion captures long-range contacts\n")
cat("  ✓ Great for TAD and boundary analysis\n")
cat("  ✗ Requires appropriate maxDist parameter\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# 6. Data Properties Summary
# ──────────────────────────────────────────────────────────────────────────────

cat("=== Data Summary ===\n\n")

cat("Contact Matrix Statistics:\n")
cat("  Total contacts:", nrow(hic_matrix), "\n")
cat("  Unique bins:", n_bins, "\n")
cat("  Contact strength (min, median, max):\n")
strength_summary <- summary(hic_matrix$strength)
cat("    Min:   ", format(strength_summary["Min."], digits = 3), "\n")
cat("    Median:", format(strength_summary["Median"], digits = 3), "\n")
cat("    Max:   ", format(strength_summary["Max."], digits = 3), "\n\n")

cat("Tiles per Style:\n")
cat("  Full:      ", nrow(hic_matrix), " tiles\n")
cat("  Diagonal:  ", sum(hic_matrix$bin_j >= hic_matrix$bin_i), " tiles\n")
cat("  Triangle:  ~", sum(hic_matrix$bin_j > hic_matrix$bin_i), " tiles (subset)\n")
cat("  Rectangle: ~", sum(hic_matrix$bin_j > hic_matrix$bin_i), " tiles (subset + expansion)\n\n")

cat("=== Visualization Complete ===\n\n")

cat("Key Features Demonstrated:\n")
cat("  ✓ Linear coordinate transformation for 45° rotation\n")
cat("  ✓ Diagonal filtering for upper triangle\n")
cat("  ✓ Window expansion for rectangle style\n")
cat("  ✓ Color mapping for contact strength\n")
cat("  ✓ SeqPlot integration with SeqTrack\n")
cat("  ✓ Genomic coordinate system handling\n\n")

cat("Next Steps:\n")
cat("  - Modify maxDist parameter to explore different scales\n")
cat("  - Try different color palettes with seq_scale_fill_continuous()\n")
cat("  - Combine multiple regions with different styles\n")
cat("  - Add annotations using SeqAnnotation elements\n")
