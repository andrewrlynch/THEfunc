# Example 03: SeqPoint with aes(color = ...) — continuous color scale
# Exercises:
#   - aes(color = vaf) + seq_scale_color_continuous(palette, limits)
#   - aes(color = vaf) with no scale (auto range)
#   - seq_scale_continuous() on the track y-axis
#
# Run with: devtools::load_all(); source("inst/examples/03_aes_continuous_color.R")

library(GenomicRanges)
devtools::load_all()

win <- CreateSequenceWindows("chr17:1-83257441")

set.seed(42)
n <- 80
snv <- GRanges(
  seqnames = "chr17",
  ranges   = IRanges(
    start = sort(sample(1:83e6, n)),
    width = 1
  ),
  vaf      = round(rbeta(n, 2, 3), 3),          # 0–1, skewed low
  depth    = sample(20:200, n, replace = TRUE)
)


# ── Continuous color on VAF ───────────────────────────────────────────────────
pt_vaf <- SeqPoint(
  snv,
  yCol  = "vaf",
  aes   = aes(color = vaf),
  scale = seq_scale_color_continuous(palette = "viridis", limits = c(0, 1))
)

# Sanity: colors are hex strings
stopifnot(all(grepl("^#[0-9A-Fa-f]{6}$", pt_vaf$aesthetics$color)))
# Low-VAF points (< 0.1) should be closer to viridis min (#440154) than max (#FDE725)
low_idx <- which(mcols(snv)$vaf < 0.1)
if (length(low_idx) > 0) {
  low_colors <- pt_vaf$aesthetics$color[low_idx]
  stopifnot(all(tolower(low_colors) != "#fde725"))
}
message("03_aes_continuous_color: viridis mapping checks passed")


# ── Plasma palette on sequencing depth ───────────────────────────────────────
pt_depth <- SeqPoint(
  snv,
  yCol  = "vaf",
  aes   = aes(color = depth),
  scale = seq_scale_color_continuous(
    palette = "plasma",
    limits  = c(20, 200),
    na_value = "grey80"
  )
)
stopifnot(all(grepl("^#[0-9A-Fa-f]{6}$", pt_depth$aesthetics$color)))
message("03_aes_continuous_color: plasma mapping checks passed")


# ── Auto-range (no explicit limits) ──────────────────────────────────────────
pt_auto <- SeqPoint(
  snv,
  yCol  = "vaf",
  aes   = aes(color = vaf)
  # no scale= → .resolve_aes() auto-ranges to range(vaf)
)
stopifnot(all(grepl("^#[0-9A-Fa-f]{6}$", pt_auto$aesthetics$color)))
message("03_aes_continuous_color: auto-range mapping checks passed")


# ── Render: two tracks, VAF colored by VAF + by depth ────────────────────────
y_scale <- seq_scale_continuous(limits = c(0, 1), n_breaks = 5)

plt <- SeqPlot(windows = win) %|%
  SeqTrack(scale_y = y_scale,
           aesthetics = list(yAxisTitleText = "VAF (viridis)")) %+% pt_vaf %|%
  SeqTrack(scale_y = y_scale,
           aesthetics = list(yAxisTitleText = "VAF (plasma = depth)")) %+% pt_depth

grid::grid.newpage()
plt$layoutGrid()
plt$drawGrid()
plt$drawAxes()
plt$drawElements()
message("03_aes_continuous_color: rendered OK")
