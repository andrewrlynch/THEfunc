# Example 07: SeqTile with style parameter — diagonal and rotated heatmaps
#
# This example demonstrates the new 'style' parameter for SeqTile, which enables
# visualization of symmetric genomic interaction matrices (Hi-C, correlation maps, etc.)
# in different formats:
#   - style="full" (default): Standard full grid of rectangles
#   - style="diagonal": Lower diagonal only (upper triangle removed)
#   - style="triangle": 45° rotated lower diagonal forming a triangle
#   - style="rectangle": 45° rotated with window expansion for symmetric appearance
#
# Key concepts:
#   - The style parameter only applies in 2D genomic mode (when y= is a GRanges)
#   - For rotated styles, tiles are drawn as diamonds using grid.polygon()
#   - The maxDist parameter controls the maximum genomic distance for clipping
#   - The yCoordType parameter specifies y-axis coordinate system ("genomic" or "distance")
#
# Run with: devtools::load_all(); source("inst/examples/07_seqtile_style_diagonal.R")

library(GenomicRanges)
library(IRanges)
library(grid)

devtools::load_all()

cat("\n=== SeqTile Style Parameter Examples ===\n")

# ══════════════════════════════════════════════════════════════════════════════
# Create sample 2D genomic data (symmetric Hi-C-style contact matrix)
# ══════════════════════════════════════════════════════════════════════════════

set.seed(42)

# Create a symmetric contact matrix for chr1:28950kb-29800kb
# This mimics a typical Hi-C heatmap region
n_bins <- 10
bin_size <- 85000  # 85 kb bins

x_starts <- seq(28950000, 28950000 + (n_bins - 1) * bin_size, by = bin_size)
y_starts <- seq(28950000, 28950000 + (n_bins - 1) * bin_size, by = bin_size)

# Create contact matrix (symmetric)
contacts <- list()
for (i in 1:n_bins) {
  for (j in i:n_bins) {
    # Contact strength decreases with distance (typical Hi-C pattern)
    distance <- abs(i - j)
    strength <- exp(-distance / 3) * runif(1, 0.5, 1.5)

    # Add contact at (i, j)
    contacts[[length(contacts) + 1]] <- list(x_idx = i, y_idx = j, strength = strength)

    # Add symmetric contact at (j, i) if i != j
    if (i != j) {
      contacts[[length(contacts) + 1]] <- list(x_idx = j, y_idx = i, strength = strength)
    }
  }
}

# Convert to GRanges
contact_df <- data.frame(
  x_idx = sapply(contacts, function(x) x$x_idx),
  y_idx = sapply(contacts, function(x) x$y_idx),
  strength = sapply(contacts, function(x) x$strength)
)

# Create x and y GRanges from contact indices
x_gr <- GRanges(
  seqnames = "chr21",
  ranges = IRanges(
    start = x_starts[contact_df$x_idx],
    width = bin_size
  ),
  strength = contact_df$strength
)

y_gr <- GRanges(
  seqnames = "chr21",
  ranges = IRanges(
    start = y_starts[contact_df$y_idx],
    width = bin_size
  ),
  strength = contact_df$strength
)

# Color by contact strength (blue = weak, red = strong)
color_palette <- colorRampPalette(c("#E8F4FF", "#4385BE", "#205EA6", "#D5300B"))(256)
colors <- color_palette[pmax(1, round(contact_df$strength * 255) + 1)]
mcols(x_gr)$color <- colors
mcols(y_gr)$color <- colors

# Define viewing window
window <- GRanges("chr21", IRanges(28950000, 29800000))
mcols(window)$scale <- 1e-6  # Display in mb

cat("Created symmetric contact matrix for chr21:28950kb-29800kb\n")
cat("  X-axis: genomic positions\n")
cat("  Y-axis: genomic positions\n")
cat("  Data: ", nrow(contact_df), " contact tiles\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Example 1: style="full" (default) — full rectangular heatmap
# ══════════════════════════════════════════════════════════════════════════════

cat("1. Creating tile with style='full' (default)...\n")

tile_full <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "full",  # Default behavior
  aesthetics = list(border = NA, lwd = 0.05)
)

stopifnot(tile_full$style == "full")
cat("   ✓ tile_full$style = ", tile_full$style, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Example 2: style="diagonal" — lower diagonal only
# ══════════════════════════════════════════════════════════════════════════════

cat("2. Creating tile with style='diagonal'...\n")

tile_diagonal <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "diagonal",
  aesthetics = list(border = NA, lwd = 0.05)
)

stopifnot(tile_diagonal$style == "diagonal")
cat("   ✓ tile_diagonal$style = ", tile_diagonal$style, "\n")
cat("   Note: diagonal filtering removes upper triangle\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Example 3: style="triangle" — rotated 45° lower diagonal forming triangle
# ══════════════════════════════════════════════════════════════════════════════

cat("3. Creating tile with style='triangle'...\n")

tile_triangle <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "triangle",
  aesthetics = list(border = NA, lwd = 0.05)
)

stopifnot(tile_triangle$style == "triangle")
stopifnot(tile_triangle$yCoordType == "distance")
cat("   ✓ tile_triangle$style = ", tile_triangle$style, "\n")
cat("   ✓ tile_triangle$yCoordType = ", tile_triangle$yCoordType, "\n")
cat("   Note: 45° rotated with tiles as diamonds\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Example 4: style="rectangle" — rotated with window expansion
# ══════════════════════════════════════════════════════════════════════════════

cat("4. Creating tile with style='rectangle'...\n")

# Use maxDist parameter to control expansion
maxDist_val <- 300000  # 300 kb expansion on each side

tile_rectangle <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "rectangle",
  maxDist = maxDist_val,
  aesthetics = list(border = NA, lwd = 0.05)
)

stopifnot(tile_rectangle$style == "rectangle")
stopifnot(tile_rectangle$maxDist == maxDist_val)
stopifnot(tile_rectangle$yCoordType == "distance")
cat("   ✓ tile_rectangle$style = ", tile_rectangle$style, "\n")
cat("   ✓ tile_rectangle$maxDist = ", tile_rectangle$maxDist, " bp\n")
cat("   ✓ tile_rectangle$yCoordType = ", tile_rectangle$yCoordType, "\n")
cat("   Note: window expanded by maxDist on both sides\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Example 5: Auto-computed maxDist for rectangle style
# ══════════════════════════════════════════════════════════════════════════════

cat("5. Creating tile with style='rectangle' and auto-computed maxDist...\n")

tile_rectangle_auto <- SeqTile(
  x = x_gr,
  y = y_gr,
  style = "rectangle"
  # maxDist defaults to max(width(y_gr))
)

cat("   ✓ Auto-computed maxDist = ", tile_rectangle_auto$maxDist, " bp\n")
cat("     (calculated as max(width(y_gr)))\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Validation tests
# ══════════════════════════════════════════════════════════════════════════════

cat("=== Validation Tests ===\n\n")

# Test 1: Invalid style reverts to full in 1D mode
cat("Test 1: Style validation\n")
tile_1d <- GRanges(
  seqnames = "chr21",
  ranges = IRanges(start = 28950000, width = 850000),
  color = "#FF0000"
)

suppressWarnings({
  tile_1d_triangle <- SeqTile(x = tile_1d, style = "triangle")
})
stopifnot(tile_1d_triangle$style == "full")
cat("  ✓ style='triangle' without y= GRanges reverts to 'full'\n\n")

# Test 2: Helper function availability
cat("Test 2: Helper functions\n")
stopifnot(exists(".validate_style_params", mode = "function"))
cat("  ✓ .validate_style_params() defined\n")
stopifnot(exists(".filter_diagonal_tiles", mode = "function"))
cat("  ✓ .filter_diagonal_tiles() defined\n")
stopifnot(exists(".transform_to_rotated_coords", mode = "function"))
cat("  ✓ .transform_to_rotated_coords() defined\n\n")

# Test 3: Linear coordinate transformation formula
cat("Test 3: Linear coordinate transformation\n")
rot <- .transform_to_rotated_coords(0, 1, 0, 1)
stopifnot(length(rot$x) == 4 && length(rot$y) == 4)
stopifnot(all(is.finite(rot$x)) && all(is.finite(rot$y)))
cat("  ✓ Unit square transformed to 4 valid corner points (diamond shape)\n")
cat("    x-coords: ", paste(round(rot$x, 3), collapse = ", "), "\n")
cat("    y-coords: ", paste(round(rot$y, 3), collapse = ", "), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Visualization Tests: Render all 4 styles with SeqPlot + SeqTrack
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Visualization Tests ===\n\n")

# Helper function to render a single plot
render_style_plot <- function(tile, style_name, window) {
  cat("Rendering style='", style_name, "'...\n")

  plt <- SeqPlot(windows = window) %|%
    SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
    tile

  grid::grid.newpage()
  plt$layoutGrid()
  plt$drawGrid()
  plt$drawAxes()
  plt$drawElements()

  invisible(plt)
}

plt <- SeqPlot(windows = window) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
  tile_full


plt$plot()

# Render each style
cat("Style 1: full (default rectangular heatmap)\n")
plt_full <- render_style_plot(tile_full, "full", window)

cat("\nStyle 2: diagonal (lower diagonal only)\n")
plt_diag <- render_style_plot(tile_diagonal, "diagonal", window)

cat("\nStyle 3: triangle (45° rotated forming triangle)\n")
plt_tri <- render_style_plot(tile_triangle, "triangle", window)

cat("\nStyle 4: rectangle (45° rotated with window expansion)\n")
plt_rect <- render_style_plot(tile_rectangle, "rectangle", window)

cat("\n=== Visualization Complete ===\n")
cat("All 4 styles rendered successfully!\n")
cat("✓ Colors represent contact strength (blue=weak, red=strong)\n")
cat("✓ Axes show genomic coordinates\n")
cat("✓ Rotated styles use diamond polygons instead of rectangles\n")
