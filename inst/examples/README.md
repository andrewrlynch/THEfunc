# SeqTile Examples

This directory contains comprehensive examples demonstrating the SeqTile element and its features, particularly the new `style` parameter for rotated heatmap visualization.

## Quick Navigation

### 📚 **Learning Path** (Recommended order)

1. **`01_operator_pipeline.R`** - Basic pipeline usage
2. **`02_aes_discrete_tile.R`** - Discrete aesthetic mapping with grouping
3. **`03_aes_continuous_color.R`** - Continuous color scales
4. **`04_position_scales.R`** - Custom position scales
5. **`05_full_multitrack.R`** - Multi-track visualization
6. **`06_multiwindow_2d_heatmap.R`** - 2D heatmaps (foundation for rotated styles)
7. **`09_style_parameter_quickstart.R`** ⭐ **START HERE** - Quick intro to new styles
8. **`07_seqtile_style_diagonal.R`** - Technical deep dive
9. **`08_rotated_hic_heatmap.R`** - Practical Hi-C visualization

## Example Descriptions

### Core Examples (Fundamentals)

#### 01_operator_pipeline.R
- **What**: Basic SeqTile usage with SeqPlot and SeqTrack
- **Topics**: Pipeline operators (`%|%`, `%+%`), window creation
- **Use case**: Getting started with SeqTile

#### 02_aes_discrete_tile.R
- **What**: Grouping tiles by discrete categories
- **Topics**: `aes(y = column)` for grouping, discrete color scales, `seq_scale_fill_discrete()`
- **Use case**: Cell type or sample visualization

#### 03_aes_continuous_color.R
- **What**: Color mapping with continuous values
- **Topics**: `fill` aesthetic, continuous scales, color ramps
- **Use case**: Expression levels, read coverage, confidence scores

#### 04_position_scales.R
- **What**: Customizing coordinate systems
- **Topics**: Custom scales, coordinate transformations, axis customization
- **Use case**: Non-linear scaling, zoomed views

#### 05_full_multitrack.R
- **What**: Combining multiple tracks with different elements
- **Topics**: Multi-element plots, track stacking, visual coordination
- **Use case**: Integrated genomic visualization

#### 06_multiwindow_2d_heatmap.R
- **What**: 2D heatmaps with genomic coordinates on both axes
- **Topics**: 2D mode with `y=GRanges`, symmetric data visualization
- **Use case**: Foundation for understanding rotated heatmaps

### Style Parameter Examples (New Features)

#### ⭐ 09_style_parameter_quickstart.R
**RECOMMENDED STARTING POINT FOR NEW USERS**

- **What**: Quick reference guide for the new `style` parameter
- **Duration**: ~15 minutes to run
- **Topics**:
  - The 4 style options (`full`, `diagonal`, `triangle`, `rectangle`)
  - When to use each style
  - Key parameters explained
  - Important implementation notes
  - Common use cases

- **Example code**:
  ```r
  # Default (full heatmap)
  tile <- SeqTile(x = x_gr, y = y_gr, style = "full")

  # Rotated triangle (most compact, standard for Hi-C)
  tile <- SeqTile(x = x_gr, y = y_gr, style = "triangle")

  # Rotated rectangle (shows off-diagonal region)
  tile <- SeqTile(x = x_gr, y = y_gr, style = "rectangle", maxDist = 2000)
  ```

#### 07_seqtile_style_diagonal.R
**TECHNICAL DEEP DIVE**

- **What**: Comprehensive technical demonstration of all style features
- **Duration**: ~30 minutes to run
- **Topics**:
  - Creating symmetric contact matrices
  - All 4 style modes in detail
  - Auto-computation of `maxDist`
  - Style validation and parameter handling
  - Helper function demonstrations
  - Rotation mathematics validation

- **Key sections**:
  - Example 1-4: Creating each style variant
  - Example 5: Auto-computed maxDist
  - Validation tests: Parameter validation, helper functions
  - Visualization tests: Full SeqPlot/SeqTrack rendering

#### 08_rotated_hic_heatmap.R
**PRACTICAL REAL-WORLD USE CASE**

- **What**: Complete Hi-C heatmap visualization workflow
- **Duration**: ~45 minutes to run
- **Topics**:
  - Realistic Hi-C contact data generation
  - Distance-dependent contact strength decay
  - Color mapping (blue=weak, red=strong)
  - Comparing all 4 visualization styles
  - TAD (Topologically Associating Domain) analysis
  - Guidance on choosing the right style

- **Example data**:
  - 100×100 bin matrix (100 kb bins)
  - chr5:40-50 Mb region
  - 10,000 contact pairs
  - Realistic decay patterns

- **Visualizations**:
  - Full heatmap (10,000 visible tiles)
  - Diagonal view (5,050 tiles, removes redundancy)
  - Triangle view (5,000 tiles, 45° rotated, compact)
  - Rectangle view (with ±2 Mb window expansion)

## The 4 Style Options

### 1. `style = "full"` (Default)

```r
tile <- SeqTile(x = x_gr, y = y_gr, style = "full")
```

**Characteristics:**
- Traditional rectangular heatmap
- Shows complete symmetric matrix
- Uses `grid.rect()` for rendering

**Advantages:**
- ✅ Shows complete picture
- ✅ Intuitive for most users
- ✅ Good for presentations

**Disadvantages:**
- ❌ Wastes space on redundant upper triangle
- ❌ Requires more vertical space

**Best for:** Publications, presentations, complete data view

---

### 2. `style = "diagonal"`

```r
tile <- SeqTile(x = x_gr, y = y_gr, style = "diagonal")
```

**Characteristics:**
- Lower diagonal only (upper triangle filtered out)
- Still uses rectangular tiles
- Uses `grid.rect()` for rendering

**Advantages:**
- ✅ Removes redundant upper triangle
- ✅ More compact than full
- ✅ Still uses rectangular tiles (intuitive)

**Disadvantages:**
- ❌ Still uses considerable space
- ❌ Not as compact as rotated styles

**Best for:** Removing redundancy while keeping familiar rectangular appearance

---

### 3. `style = "triangle"` (45° Rotated)

```r
tile <- SeqTile(x = x_gr, y = y_gr, style = "triangle")
```

**Characteristics:**
- Linear coordinate transformation: `x_rot = (x+y)/2`, `y_rot = (y-x)/2`
- Lower diagonal rotated 45° forming triangle
- Uses `grid.polygon()` for rendering diamonds
- Most compact representation

**Advantages:**
- ✅ Most space-efficient
- ✅ Standard in bioinformatics (Hi-C visualization)
- ✅ Shows local TAD structure clearly
- ✅ Professional appearance

**Disadvantages:**
- ❌ Less intuitive for unfamiliar viewers
- ❌ Requires explanation in figures

**Best for:** Hi-C heatmaps, publication-quality figures, space-constrained layouts

---

### 4. `style = "rectangle"` (45° Rotated with Expansion)

```r
tile <- SeqTile(x = x_gr, y = y_gr, style = "rectangle", maxDist = 2000)
```

**Characteristics:**
- Linear coordinate transformation (same as triangle)
- Rotated lower diagonal in rectangular frame
- Window expansion shows off-diagonal contacts
- Configurable expansion distance (`maxDist`)

**Parameters:**
- `maxDist = NULL` (default): Auto-computed from max y-range width
- `maxDist = <numeric>`: Custom expansion distance in bp

**Advantages:**
- ✅ Shows off-diagonal region
- ✅ Window expansion captures long-range contacts
- ✅ Great for TAD boundary analysis
- ✅ Configurable to explore different scales

**Disadvantages:**
- ❌ Requires appropriate `maxDist` selection
- ❌ May show empty space if maxDist too large

**Best for:** TAD analysis, long-range interaction studies, customizable visualization scales

---

## Key Parameters Reference

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `x` | GRanges | required | Genomic ranges for x-axis |
| `y` | GRanges | NULL | Genomic ranges for y-axis (2D mode) |
| `style` | character | "full" | Visualization style: "full", "diagonal", "triangle", "rectangle" |
| `maxDist` | numeric | NULL | Maximum distance for rectangle style (auto-computed if NULL) |

### Other Parameters

- `aes`: Aesthetic mappings (`y` for grouping, `fill` for color)
- `scale`: Color scale (`seq_scale_fill_discrete()`, etc.)
- `aesthetics`: Visual parameters (`border`, `lwd`)

## Running the Examples

### Run a single example:
```r
devtools::load_all()
source("inst/examples/09_style_parameter_quickstart.R")
```

### Run all examples:
```bash
for i in 01 02 03 04 05 06 07 08 09; do
  echo "Running example $i..."
  R --slave -e "devtools::load_all(); source('inst/examples/${i:0:2}_*.R')"
done
```

### Run in RStudio:
1. Open RStudio
2. Go to File → Open File
3. Select the example file
4. Press Ctrl+Shift+S (or Cmd+Shift+S on Mac) to source the entire file
5. View plots as they render

## Important Notes

### 2D Mode Requirement
- Style parameter **only works when `y=` is provided as GRanges**
- 1D mode (single axis): style always reverts to "full"
- Useful error message warns if style is used incorrectly

### Mathematical Foundation
- Rotated styles use **linear coordinate transformation**, not geometric rotation
- Transformation: `x_rot = (x + y) / 2`, `y_rot = (y - x) / 2`
- This matches standard Hi-C visualization approach
- Based on proven implementation from `hic_rotated_heatmap.R`

### Backward Compatibility
- Default `style = "full"` preserves existing behavior exactly
- All existing code works unchanged
- No breaking changes to SeqTile API

### Performance Notes
- `style = "full"`: Fast, uses `grid.rect()` vectorization
- `style = "diagonal"`: Very fast, minimal filtering overhead
- `style = "triangle"`: Fast, uses vectorized filtering + rotation
- `style = "rectangle"`: Fast, includes window expansion (negligible overhead)

## Troubleshooting

### Plot appears blank
- Check that `x=` is provided as GRanges
- Ensure genomic ranges overlap with viewing window
- Verify window is created with `createGenomeWindows()` or `GRanges()`

### Tiles don't appear in rotated styles
- Verify `y=` is provided (required for rotated styles)
- Check that `maxDist` is reasonable for your data scale
- Ensure yCoordType matches your axis configuration

### Style parameter not recognized
- Make sure you have latest version: `devtools::load_all()`
- Verify `y=` parameter is provided (style requires 2D mode)
- Check package loaded correctly: `library(THEfunc)`

## Further Reading

- **SeqPlot**: See examples 01, 05 for SeqPlot usage
- **SeqTrack**: See examples 02, 06 for track configuration
- **Color scales**: See examples 03, 04 for aesthetic and scale options
- **Genomic data**: Read Bioconductor's GRanges documentation

## Citation

If you use these examples in your work, please cite:
- THEfunc package
- Bioconductor GenomicRanges package

## Questions?

Refer to the specific example's comments for detailed explanations:
- Search for `#' @description` for method documentation
- Look for `cat()` statements for narrative explanations
- Check helper function definitions for implementation details
