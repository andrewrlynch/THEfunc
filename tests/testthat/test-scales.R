# ── Color scales ──────────────────────────────────────────────────────────────

test_that("seq_scale_color_continuous() returns correct class", {
  s <- seq_scale_color_continuous()
  expect_s3_class(s, "SeqScaleContinuous")
  expect_s3_class(s, "SeqScale")
})

test_that("seq_scale_color_continuous() defaults are correct", {
  s <- seq_scale_color_continuous()
  expect_equal(s$aesthetic, "color")
  expect_equal(s$type, "continuous")
  expect_equal(s$palette, "viridis")
  expect_null(s$limits)
  expect_equal(s$na_value, "grey80")
})

test_that("seq_scale_color_continuous() accepts custom palette and limits", {
  s <- seq_scale_color_continuous(palette = "plasma", limits = c(0, 100), na_value = "#FF0000")
  expect_equal(s$palette, "plasma")
  expect_equal(s$limits, c(0, 100))
  expect_equal(s$na_value, "#FF0000")
})

test_that("seq_scale_color_discrete() returns correct class", {
  s <- seq_scale_color_discrete()
  expect_s3_class(s, "SeqScaleDiscrete")
  expect_s3_class(s, "SeqScale")
})

test_that("seq_scale_color_discrete() defaults are correct", {
  s <- seq_scale_color_discrete()
  expect_equal(s$aesthetic, "color")
  expect_equal(s$type, "discrete")
  expect_null(s$values)
  expect_null(s$palette)
  expect_equal(s$na_value, "grey80")
})

test_that("seq_scale_color_discrete() stores named color values", {
  vals <- c(TypeA = "#FF0000", TypeB = "#0000FF")
  s <- seq_scale_color_discrete(values = vals)
  expect_equal(s$values, vals)
})

test_that("seq_scale_color_discrete() stores a palette function", {
  pal_fn <- function(n) grDevices::rainbow(n)
  s <- seq_scale_color_discrete(palette = pal_fn)
  expect_identical(s$palette, pal_fn)
})

# ── Fill scale delegation ─────────────────────────────────────────────────────

test_that("seq_scale_fill_continuous() sets aesthetic to 'fill'", {
  s <- seq_scale_fill_continuous()
  expect_equal(s$aesthetic, "fill")
  expect_s3_class(s, "SeqScaleContinuous")
})

test_that("seq_scale_fill_continuous() passes through all other args", {
  s <- seq_scale_fill_continuous(palette = "magma", limits = c(-1, 1))
  expect_equal(s$palette, "magma")
  expect_equal(s$limits, c(-1, 1))
})

test_that("seq_scale_fill_discrete() sets aesthetic to 'fill'", {
  s <- seq_scale_fill_discrete()
  expect_equal(s$aesthetic, "fill")
  expect_s3_class(s, "SeqScaleDiscrete")
})

test_that("seq_scale_fill_discrete() passes through values", {
  vals <- c(A = "red", B = "blue")
  s <- seq_scale_fill_discrete(values = vals)
  expect_equal(s$values, vals)
})

# ── Position scales ───────────────────────────────────────────────────────────

test_that("seq_scale_genomic() requires a GRanges object", {
  expect_error(seq_scale_genomic("not_granges"))
  expect_error(seq_scale_genomic(data.frame(chr = "chr1")))
})

test_that("seq_scale_genomic() returns correct class", {
  gr <- GenomicRanges::GRanges("chr1:1-1000000")
  s <- seq_scale_genomic(gr)
  expect_s3_class(s, "SeqScaleGenomic")
  expect_s3_class(s, "SeqPositionScale")
})

test_that("seq_scale_genomic() stores windows and defaults scale_factor to 1e-6", {
  gr <- GenomicRanges::GRanges("chr1:1-1000000")
  s <- seq_scale_genomic(gr)
  expect_identical(s$windows, gr)
  expect_equal(s$type, "genomic")
  expect_equal(s$scale_factor, 1e-6)
})

test_that("seq_scale_genomic() reads scale_factor from mcols$scale", {
  gr <- GenomicRanges::GRanges("chr1:1-1000000")
  S4Vectors::mcols(gr)$scale <- 1e-3
  s <- seq_scale_genomic(gr)
  expect_equal(s$scale_factor, 1e-3)
})

test_that("seq_scale_genomic() accepts explicit scale_factor override", {
  gr <- GenomicRanges::GRanges("chr1:1-1000000")
  s <- seq_scale_genomic(gr, scale_factor = 1)
  expect_equal(s$scale_factor, 1)
})

test_that("seq_scale_continuous() returns correct class and defaults", {
  s <- seq_scale_continuous()
  expect_s3_class(s, "SeqScaleContinuous_Pos")
  expect_s3_class(s, "SeqPositionScale")
  expect_equal(s$type, "continuous")
  expect_null(s$limits)
  expect_equal(s$n_breaks, 5)
})

test_that("seq_scale_continuous() accepts limits and n_breaks", {
  s <- seq_scale_continuous(limits = c(0, 1), n_breaks = 10)
  expect_equal(s$limits, c(0, 1))
  expect_equal(s$n_breaks, 10)
})

test_that("seq_scale_discrete() returns correct class and defaults", {
  s <- seq_scale_discrete()
  expect_s3_class(s, "SeqScaleDiscrete_Pos")
  expect_s3_class(s, "SeqPositionScale")
  expect_equal(s$type, "discrete")
  expect_null(s$levels)
  expect_null(s$labels)
})

test_that("seq_scale_discrete() stores levels and labels", {
  s <- seq_scale_discrete(levels = c("T_cell", "B_cell"), labels = c("T", "B"))
  expect_equal(s$levels, c("T_cell", "B_cell"))
  expect_equal(s$labels, c("T", "B"))
})
