# Example 09: SeqTile Style Parameter - Quick Start Guide
#
# A quick reference guide for using the new style parameter in SeqTile.
# Perfect for getting started with rotated heatmap visualization.
#
# Run with: devtools::load_all(); source("inst/examples/09_style_parameter_quickstart.R")

library(GenomicRanges)
devtools::load_all()

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║      SeqTile Style Parameter - Quick Start Guide          ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Quick Setup: Create Sample 2D Data
# ──────────────────────────────────────────────────────────────────────────────

cat("STEP 1: Create Sample Data\n")
cat("────────────────────────────────────────────────────────────\n\n")

set.seed(42)
n <- 6  # 6x6 matrix of contacts

# Create all contact pairs for 6x6 matrix
contact_df <- expand.grid(x_idx = 1:n, y_idx = 1:n)
bin_starts <- seq(1000, 6000, by = 1000)
bin_size <- 500

# Create strength values (decay with distance)
strength <- exp(-abs(contact_df$x_idx - contact_df$y_idx) / 2) * runif(nrow(contact_df), 0.5, 1.5)

# Create color palette and assign colors based on strength
color_palette <- colorRampPalette(c("#E8F4FF", "#4385BE", "#205EA6", "#D5300B"))(256)
colors <- color_palette[pmax(1, round(strength / max(strength) * 255) + 1)]

# Create genomic ranges for x-axis (each contact pair references a bin)
x_gr <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = bin_starts[contact_df$x_idx], width = bin_size),
  color = colors
)

# Create genomic ranges for y-axis (each contact pair references a bin)
y_gr <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = bin_starts[contact_df$y_idx], width = bin_size),
  color = colors
)

# Create viewing window
window <- GRanges("chr1", IRanges(1000, 6500))

cat("✓ Created", n, "×", n, "contact matrix\n")
cat("✓ Genomic region: chr1:1000-6500\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Using Different Styles
# ──────────────────────────────────────────────────────────────────────────────

cat("\nSTEP 2: Using Different Styles\n")
cat("────────────────────────────────────────────────────────────\n\n")

# Style 1: Default (full heatmap)
cat("1️⃣  FULL STYLE (Default)\n")
cat("   SeqTile(x = x_gr, y = y_gr)  # or style = \"full\"\n")
cat("   → Shows all contacts in rectangular grid\n\n")

tile_full <- SeqTile(x = x_gr, y = y_gr, style = "full")
cat("   Created: tile_full$style = \"", tile_full$style, "\"\n\n")

# Style 2: Diagonal (upper triangle only)
cat("2️⃣  DIAGONAL STYLE\n")
cat("   SeqTile(x = x_gr, y = y_gr, style = \"diagonal\")\n")
cat("   → Shows only upper triangle (removes redundancy)\n\n")

tile_diag <- SeqTile(x = x_gr, y = y_gr, style = "diagonal")
cat("   Created: tile_diag$style = \"", tile_diag$style, "\"\n\n")

# Style 3: Triangle (45° rotated, most compact)
cat("3️⃣  TRIANGLE STYLE\n")
cat("   SeqTile(x = x_gr, y = y_gr, style = \"triangle\")\n")
cat("   → 45° rotated lower diagonal forming triangle\n")
cat("   → Standard for Hi-C visualization\n\n")

tile_tri <- SeqTile(x = x_gr, y = y_gr, style = "triangle")
cat("   Created: tile_tri$style = \"", tile_tri$style, "\"\n")
cat("   yCoordType = \"", tile_tri$yCoordType, "\" (distance coordinates)\n\n")

# Style 4: Rectangle (rotated with window expansion)
cat("4️⃣  RECTANGLE STYLE\n")
cat("   SeqTile(x = x_gr, y = y_gr, style = \"rectangle\", maxDist = 1000)\n")
cat("   → 45° rotated in rectangular frame with window expansion\n")
cat("   → Shows off-diagonal contacts from expanded region\n\n")

tile_rect <- SeqTile(x = x_gr, y = y_gr, style = "rectangle", maxDist = 1000)
cat("   Created: tile_rect$style = \"", tile_rect$style, "\"\n")
cat("   maxDist = ", tile_rect$maxDist, " bp (±1000 bp expansion)\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Rendering with SeqPlot
# ──────────────────────────────────────────────────────────────────────────────

cat("\nSTEP 3: Rendering with SeqPlot\n")
cat("────────────────────────────────────────────────────────────\n\n")

cat("Creating plot with style = \"full\"...\n")

# Create and render plot
plt <- SeqPlot(windows = window) %|%
  SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
  tile_full

grid::grid.newpage()
plt$layoutGrid()
plt$drawGrid()
plt$drawAxes()
plt$drawElements()

cat("✓ Plot rendered\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Key Parameters Explained
# ──────────────────────────────────────────────────────────────────────────────

cat("\nKEY PARAMETERS\n")
cat("────────────────────────────────────────────────────────────\n\n")

cat("style:\n")
cat("  • \"full\"      → Full rectangular heatmap [DEFAULT]\n")
cat("  • \"diagonal\"  → Upper diagonal rectangle\n")
cat("  • \"triangle\"  → 45° rotated forming triangle\n")
cat("  • \"rectangle\" → 45° rotated in rectangular frame\n\n")

cat("yCoordType:  (Optional, stored for display interpretation)\n")
cat("  • \"genomic\"   → Y-axis shows genomic coordinates [DEFAULT]\n")
cat("  • \"distance\"  → Y-axis shows genomic distance\n")
cat("  Note: Only relevant for rotated styles\n\n")

cat("maxDist:  (Optional, for rectangle style)\n")
cat("  • NULL          → Auto-computed from y-axis window widths [DEFAULT]\n")
cat("  • Numeric value → Set custom expansion distance (in bp)\n")
cat("  Note: Controls how much window is expanded for rectangle style\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Important Notes
# ──────────────────────────────────────────────────────────────────────────────

cat("\nIMPORTANT NOTES\n")
cat("────────────────────────────────────────────────────────────\n\n")

cat("✓ style parameter ONLY works in 2D genomic mode:\n")
cat("  • Must provide BOTH x= and y= as GRanges\n")
cat("  • If y= is NULL, style reverts to \"full\"\n\n")

cat("✓ Rotated styles (triangle, rectangle) use linear coordinate\n")
cat("  transformation, not geometric rotation:\n")
cat("  • x_rot = (x + y) / 2\n")
cat("  • y_rot = (y - x) / 2\n\n")

cat("✓ Backward compatible:\n")
cat("  • Default style=\"full\" preserves existing behavior\n")
cat("  • All existing code works unchanged\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Common Use Cases
# ──────────────────────────────────────────────────────────────────────────────

cat("\nCOMMON USE CASES\n")
cat("────────────────────────────────────────────────────────────\n\n")

cat("Hi-C Contact Maps:\n")
cat("  → Use style=\"triangle\" for publication-quality figures\n")
cat("  → Use style=\"rectangle\" to examine TAD boundaries\n\n")

cat("Correlation Matrices:\n")
cat("  → Use style=\"full\" for complete view\n")
cat("  → Use style=\"diagonal\" if size is constraining\n\n")

cat("Genome Interaction Networks:\n")
cat("  → Use style=\"rectangle\" with custom maxDist\n")
cat("  → Adjust maxDist to explore different interaction scales\n\n")

# ──────────────────────────────────────────────────────────────────────────────
# Visual Comparison
# ──────────────────────────────────────────────────────────────────────────────

cat("\nVISUAL COMPARISON\n")
cat("────────────────────────────────────────────────────────────\n\n")

cat("Rendering all 4 styles for comparison...\n\n")

styles <- list(
  list(name = "full", tile = tile_full, desc = "Full rectangular heatmap"),
  list(name = "diagonal", tile = tile_diag, desc = "Upper diagonal only"),
  list(name = "triangle", tile = tile_tri, desc = "45° rotated triangle"),
  list(name = "rectangle", tile = tile_rect, desc = "45° rotated rectangle")
)

for (i in seq_along(styles)) {
  cat(i, ".", styles[[i]]$name, ": ", styles[[i]]$desc, "\n")

  plt <- SeqPlot(windows = window) %|%
    SeqTrack(aesthetics = list(windowBorder = NA, windowBackground = NA)) %+%
    styles[[i]]$tile

  grid::grid.newpage()
  plt$layoutGrid()
  plt$drawGrid()
  plt$drawAxes()
  plt$drawElements()
}

# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║                     QUICK START COMPLETE                  ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("You've learned:\n")
cat("  ✓ How to use the style parameter in SeqTile\n")
cat("  ✓ The 4 different visualization styles\n")
cat("  ✓ Key parameters: style, yCoordType, maxDist\n")
cat("  ✓ How to render plots with SeqPlot + SeqTrack\n\n")

cat("Next steps:\n")
cat("  1. Read example 08_rotated_hic_heatmap.R for realistic use case\n")
cat("  2. Read example 07_seqtile_style_diagonal.R for detailed technical info\n")
cat("  3. Experiment with different maxDist values\n")
cat("  4. Combine with color scales for enhanced visualization\n\n")

cat("Happy plotting! 🎨📊\n\n")
