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
