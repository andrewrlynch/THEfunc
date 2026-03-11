# Helpers ─────────────────────────────────────────────────────────────────────

# Three ranges spanning chr1, annotated with cell_type, mutation_count, and
# a pre-existing 'color' column for legacy-path tests.
.tile_gr <- function() {
  GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200", "chr1:201-300"),
    cell_type      = c("T_cell", "B_cell", "T_cell"),
    mutation_count = c(0L, 50L, 100L),
    color          = c("#FF0000", "#00FF00", "#FF0000")
  )
}

.point_gr <- function() {
  GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200", "chr1:201-300"),
    score = c(1.0, 2.0, 3.0),
    group = c("A", "B", "A")
  )
}


# ── Constructor wrappers ──────────────────────────────────────────────────────

test_that("SeqPlot() wrapper creates an R6 SeqPlot object", {
  plt <- SeqPlot()
  expect_true(inherits(plt, "SeqPlot"))
  expect_true(R6::is.R6(plt))
})

test_that("SeqTrack() wrapper creates an R6 SeqTrack object", {
  trk <- SeqTrack()
  expect_true(inherits(trk, "SeqTrack"))
  expect_true(R6::is.R6(trk))
})

test_that("SeqPoint() wrapper creates an R6 SeqPoint object", {
  gr  <- .point_gr()
  pt  <- SeqPoint(gr, yCol = "score")
  expect_true(inherits(pt, "SeqPoint"))
  expect_true(R6::is.R6(pt))
})

test_that("SeqLine() wrapper creates an R6 SeqLine object", {
  gr  <- .point_gr()
  ln  <- SeqLine(gr, yCol = "score")
  expect_true(inherits(ln, "SeqLine"))
})

test_that("SeqBar() wrapper creates an R6 SeqBar object", {
  gr  <- .point_gr()
  br  <- SeqBar(gr, yCol = "score")
  expect_true(inherits(br, "SeqBar"))
})


# ── SeqTile: new aes(y = ...) path ───────────────────────────────────────────

test_that("SeqTile aes(y=) stores groups in order of first appearance", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type))
  expect_equal(tile$groups, c("T_cell", "B_cell"))
})

test_that("SeqTile aes(y=) computes y index mapping correctly", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type))
  # T_cell=1, B_cell=2, T_cell=1
  expect_equal(tile$y, c(1L, 2L, 1L))
})

test_that("SeqTile aes(y=) sets groupCol to the mapped column name", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type))
  expect_equal(tile$groupCol, "cell_type")
})

test_that("SeqTile aes(y=) with a factor preserves factor levels", {
  gr <- .tile_gr()
  GenomicRanges::mcols(gr)$cell_type <- factor(
    GenomicRanges::mcols(gr)$cell_type,
    levels = c("B_cell", "T_cell")
  )
  tile <- SeqTile(gr, aes = aes(y = cell_type))
  expect_equal(tile$groups, c("B_cell", "T_cell"))  # factor order respected
})

test_that("SeqTile aes(y=) errors when mapped column is absent", {
  gr <- .tile_gr()
  expect_error(
    SeqTile(gr, aes = aes(y = nonexistent_col)),
    "not found"
  )
})

test_that("SeqTile aes(fill=) resolves continuous fill into mcols$color", {
  scale <- seq_scale_fill_continuous(palette = "viridis",
                                     limits  = c(0, 100))
  tile  <- SeqTile(.tile_gr(),
                   aes   = aes(y = cell_type, fill = mutation_count),
                   scale = scale)
  colors <- as.character(GenomicRanges::mcols(tile$gr)$color)
  expect_length(colors, 3L)
  expect_match(colors[1], "^#[0-9A-Fa-f]{6}$")
  # Observation at count=0 and count=100 should have different colors
  expect_false(colors[1] == colors[3])
})

test_that("SeqTile aes(fill=) viridis min/max colors are correct", {
  gr    <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200"),
    cell  = c("A", "B"),
    val   = c(0L, 100L)
  )
  scale <- seq_scale_fill_continuous(palette = "viridis", limits = c(0, 100))
  tile  <- SeqTile(gr, aes = aes(y = cell, fill = val), scale = scale)
  colors <- as.character(GenomicRanges::mcols(tile$gr)$color)
  expect_equal(toupper(colors[1]), "#440154")  # viridis low
  expect_equal(toupper(colors[2]), "#FDE725")  # viridis high
})

test_that("SeqTile aes(fill=) discrete fill uses SeqScaleDiscrete colors", {
  scale <- seq_scale_fill_discrete(values = c(T_cell = "#AA0000",
                                               B_cell = "#0000AA"))
  tile  <- SeqTile(.tile_gr(),
                   aes   = aes(y = cell_type, fill = cell_type),
                   scale = scale)
  colors <- as.character(GenomicRanges::mcols(tile$gr)$color)
  expect_equal(colors[1], "#AA0000")  # T_cell
  expect_equal(colors[2], "#0000AA")  # B_cell
  expect_equal(colors[3], "#AA0000")  # T_cell again
})

test_that("SeqTile .infer_scale_y() returns a SeqScaleDiscrete_Pos when groups are set", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type))
  scale <- tile$.infer_scale_y()
  expect_s3_class(scale, "SeqScaleDiscrete_Pos")
  expect_equal(scale$levels, tile$groups)
})

test_that("SeqTile .infer_scale_y() returns NULL when no grouping is set", {
  # A tile with only a color column, no y mapping
  gr   <- .tile_gr()
  tile <- SeqTile(gr)
  expect_null(tile$.infer_scale_y())
})


# ── SeqTile: legacy groupCol path ────────────────────────────────────────────

test_that("SeqTile groupCol (legacy) sets groups and y correctly", {
  tile <- SeqTile(.tile_gr(), groupCol = "cell_type")
  expect_equal(tile$groups, c("T_cell", "B_cell"))
  expect_equal(tile$y, c(1L, 2L, 1L))
})

test_that("SeqTile groupCol (legacy) errors when column is absent", {
  gr <- .tile_gr()
  expect_error(SeqTile(gr, groupCol = "no_such_col"), "not found")
})

test_that("SeqTile without color column and without aes(fill=) errors", {
  gr <- GenomicRanges::GRanges(
    c("chr1:1-100"),
    cell_type = "T_cell"
  )
  # No 'color' mcol and no fill aes → should error
  expect_error(SeqTile(gr, groupCol = "cell_type"))
})


# ── SeqPoint: aes + scale ────────────────────────────────────────────────────

test_that("SeqPoint with aes(color=) resolves discrete colors", {
  gr  <- .point_gr()
  sc  <- seq_scale_color_discrete(values = c(A = "#FF0000", B = "#0000FF"))
  pt  <- SeqPoint(gr, yCol = "score", aes = aes(color = group), scale = sc)
  # Obs 1 = A, 2 = B, 3 = A
  expect_equal(pt$aesthetics$color[1], "#FF0000")
  expect_equal(pt$aesthetics$color[2], "#0000FF")
  expect_equal(pt$aesthetics$color[3], "#FF0000")
})

test_that("SeqPoint with aes(color=) resolves continuous colors", {
  gr  <- .point_gr()
  sc  <- seq_scale_color_continuous(palette = "viridis",
                                    limits  = c(1, 3))
  pt  <- SeqPoint(gr, yCol = "score", aes = aes(color = score), scale = sc)
  colors <- pt$aesthetics$color
  expect_length(colors, 3L)
  expect_match(colors[1], "^#[0-9A-Fa-f]{6}$")
  # score=1 (min) and score=3 (max) should differ
  expect_false(colors[1] == colors[3])
})

test_that("SeqPoint without aes uses fixed color from aesthetics", {
  gr <- .point_gr()
  pt <- SeqPoint(gr, yCol = "score", color = "#123456")
  expect_true(all(pt$aesthetics$color == "#123456"))
})


# ── SeqSegment: aes + scale ──────────────────────────────────────────────────

test_that("SeqSegment with aes(color=) resolves colors", {
  gr <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:200-300"),
    y0   = c(0.0, 0.5),
    y1   = c(0.4, 0.9),
    grp  = c("up", "down")
  )
  sc  <- seq_scale_color_discrete(values = c(up = "#00CC00", down = "#CC0000"))
  seg <- SeqSegment(gr, y0Col = "y0", y1Col = "y1",
                    aes = aes(color = grp), scale = sc)
  expect_equal(seg$aesthetics$color[1], "#00CC00")
  expect_equal(seg$aesthetics$color[2], "#CC0000")
})


# ── SeqBar: aes + scale ──────────────────────────────────────────────────────

test_that("SeqBar with aes(fill=) resolves fill", {
  gr <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200"),
    height = c(10.0, 20.0),
    grp    = c("A", "B")
  )
  sc <- seq_scale_fill_discrete(values = c(A = "#AAAA00", B = "#0000AA"))
  br <- SeqBar(gr, yCol = "height", aes = aes(fill = grp), scale = sc)
  expect_equal(br$aesthetics$fill[1], "#AAAA00")
  expect_equal(br$aesthetics$fill[2], "#0000AA")
})


# ── SeqLine: aes + scale ─────────────────────────────────────────────────────

test_that("SeqLine with aes(color=) resolves color", {
  gr <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200", "chr1:201-300"),
    val  = c(1.0, 5.0, 10.0),
    grp  = c("X", "X", "X")
  )
  sc <- seq_scale_color_continuous(palette = "plasma", limits = c(1, 10))
  ln <- SeqLine(gr, yCol = "val", aes = aes(color = val), scale = sc)
  expect_length(ln$aesthetics$color, 3L)
  expect_match(ln$aesthetics$color[1], "^#[0-9A-Fa-f]{6}$")
})


# ── SeqTile: style parameter tests ────────────────────────────────────────────

test_that("SeqTile style parameter defaults to 'full'", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type))
  expect_equal(tile$style, "full")
})

test_that("SeqTile style parameter accepts valid values", {
  gr_x <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200", "chr1:201-300"),
    color = c("#FF0000", "#00FF00", "#0000FF")
  )
  gr_y <- GenomicRanges::GRanges(
    c("chr1:50-150", "chr1:150-250", "chr1:250-350"),
    color = c("#FF0000", "#00FF00", "#0000FF")
  )

  tile_full <- SeqTile(x = gr_x, y = gr_y, style = "full")
  expect_equal(tile_full$style, "full")

  tile_diag <- SeqTile(x = gr_x, y = gr_y, style = "diagonal")
  expect_equal(tile_diag$style, "diagonal")

  tile_tri <- SeqTile(x = gr_x, y = gr_y, style = "triangle")
  expect_equal(tile_tri$style, "triangle")

  tile_rect <- SeqTile(x = gr_x, y = gr_y, style = "rectangle")
  expect_equal(tile_rect$style, "rectangle")
})

test_that("SeqTile invalid style value raises error", {
  gr_x <- GenomicRanges::GRanges(
    c("chr1:1-100"),
    color = c("#FF0000")
  )
  gr_y <- GenomicRanges::GRanges(
    c("chr1:50-150"),
    color = c("#FF0000")
  )

  expect_error(
    SeqTile(x = gr_x, y = gr_y, style = "invalid"),
    "style must be one of"
  )
})

test_that("SeqTile style != 'full' with no gr_y reverts to 'full' with warning", {
  tile <- SeqTile(.tile_gr(), aes = aes(y = cell_type), style = "triangle")
  expect_equal(tile$style, "full")
  expect_warning(
    SeqTile(.tile_gr(), aes = aes(y = cell_type), style = "triangle"),
    "Reverting to 'full'"
  )
})

test_that("SeqTile maxDist parameter accepts positive numeric values", {
  gr_x <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200"),
    color = c("#FF0000", "#00FF00")
  )
  gr_y <- GenomicRanges::GRanges(
    c("chr1:50-150", "chr1:150-250"),
    color = c("#FF0000", "#00FF00")
  )

  tile <- SeqTile(x = gr_x, y = gr_y, style = "rectangle", maxDist = 50)
  expect_equal(tile$maxDist, 50)
})

test_that("SeqTile maxDist must be positive", {
  gr_x <- GenomicRanges::GRanges(c("chr1:1-100"), color = "#FF0000")
  gr_y <- GenomicRanges::GRanges(c("chr1:50-150"), color = "#FF0000")

  expect_error(
    SeqTile(x = gr_x, y = gr_y, style = "rectangle", maxDist = -50),
    "must be a positive number"
  )
})

test_that("SeqTile rectangle style auto-computes maxDist from gr_y width", {
  gr_x <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200"),
    color = c("#FF0000", "#00FF00")
  )
  gr_y <- GenomicRanges::GRanges(
    c("chr1:50-150", "chr1:150-350"),  # widths: 100, 200
    color = c("#FF0000", "#00FF00")
  )

  tile <- SeqTile(x = gr_x, y = gr_y, style = "rectangle")
  expect_equal(tile$maxDist, 200)  # max width
})

test_that("SeqTile yCoordType parameter is stored", {
  gr_x <- GenomicRanges::GRanges(c("chr1:1-100"), color = "#FF0000")
  gr_y <- GenomicRanges::GRanges(c("chr1:50-150"), color = "#FF0000")

  tile_genomic <- SeqTile(x = gr_x, y = gr_y, yCoordType = "genomic")
  expect_equal(tile_genomic$yCoordType, "genomic")

  tile_dist <- SeqTile(x = gr_x, y = gr_y, yCoordType = "distance")
  expect_equal(tile_dist$yCoordType, "distance")
})

test_that("SeqTile rotation function .rotate_coordinates_45 works correctly", {
  # Test a unit square [0, 1] x [0, 1]
  rot <- .rotate_coordinates_45(0, 1, 0, 1)

  expect_length(rot$x, 4)
  expect_length(rot$y, 4)

  # Rotation should move corners symmetrically
  # Original corners: (0,0), (1,0), (1,1), (0,1)
  # Center: (0.5, 0.5)
  cos45 <- sqrt(2) / 2
  sin45 <- sqrt(2) / 2

  # Expected bottom-left corner rotation: (-0.5, -0.5) rotated
  expected_x1 <- 0.5 + (-0.5) * cos45 - (-0.5) * sin45
  expected_y1 <- 0.5 + (-0.5) * sin45 + (-0.5) * cos45

  expect_equal(rot$x[1], expected_x1, tolerance = 1e-10)
  expect_equal(rot$y[1], expected_y1, tolerance = 1e-10)

  # Check that rotated corners are finite
  expect_true(all(is.finite(rot$x)))
  expect_true(all(is.finite(rot$y)))
})

test_that("SeqTile coordCanvas includes original coordinates for 2D mode", {
  gr_x <- GenomicRanges::GRanges(
    c("chr1:1-100", "chr1:101-200"),
    color = c("#FF0000", "#00FF00")
  )
  gr_y <- GenomicRanges::GRanges(
    c("chr1:50-150", "chr1:150-250"),
    color = c("#FF0000", "#00FF00")
  )

  tile <- SeqTile(x = gr_x, y = gr_y, style = "diagonal")

  # After prep(), coordCanvas should have original coordinate columns
  track_windows <- gr_x
  layout_track <- list(
    list(
      xscale = c(0, 300),
      yscale = c(0, 300),
      inner = list(x0 = 0, x1 = 1, y0 = 0, y1 = 1),
      y_sub_panels = NULL
    )
  )

  tile$prep(layout_track, track_windows)

  if (length(tile$coordCanvas) > 0 && !is.null(tile$coordCanvas[[1]])) {
    coords <- tile$coordCanvas[[1]]
    expect_true("x0_orig" %in% names(coords))
    expect_true("x1_orig" %in% names(coords))
    expect_true("y0_orig" %in% names(coords))
    expect_true("y1_orig" %in% names(coords))
  }
})
