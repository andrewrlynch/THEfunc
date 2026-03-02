# Helper: a minimal single-window GRanges for SeqPlot/SeqTrack setup
.make_windows <- function() GenomicRanges::GRanges("chr1:1-248956422")

# Helper: a small GRanges with a numeric score column
.make_gr <- function() {
  GenomicRanges::GRanges(
    c("chr1:1000-2000", "chr1:3000-4000", "chr1:5000-6000"),
    score = c(1.0, 2.0, 3.0)
  )
}


# ── %|% operator ─────────────────────────────────────────────────────────────

test_that("%|% adds a SeqTrack to an empty SeqPlot", {
  plt <- SeqPlot(windows = .make_windows())
  trk <- SeqTrack()
  result <- plt %|% trk
  expect_identical(result, plt)
  expect_length(plt$tracks, 1L)
})

test_that("%|% adds multiple tracks in order", {
  plt <- SeqPlot(windows = .make_windows())
  plt %|% SeqTrack() %|% SeqTrack()
  expect_length(plt$tracks, 2L)
})

test_that("%|% returns the SeqPlot invisibly (same object)", {
  plt <- SeqPlot(windows = .make_windows())
  result <- plt %|% SeqTrack()
  expect_identical(result, plt)
})

test_that("%|% errors when right-hand side is not a SeqTrack", {
  plt <- SeqPlot(windows = .make_windows())
  expect_error(plt %|% "not a track")
  expect_error(plt %|% 42)
  expect_error(plt %|% list())
})


# ── %+% on SeqPlot ───────────────────────────────────────────────────────────

test_that("%+% adds a SeqElement to the last track in a SeqPlot", {
  plt <- SeqPlot(windows = .make_windows()) %|% SeqTrack()
  pt  <- SeqPoint(.make_gr(), yCol = "score")
  result <- plt %+% pt
  expect_identical(result, plt)
  expect_length(plt$tracks[[1]]$elements, 1L)
})

test_that("%+% targets only the last track, not earlier ones", {
  plt <- SeqPlot(windows = .make_windows()) %|% SeqTrack() %|% SeqTrack()
  pt  <- SeqPoint(.make_gr(), yCol = "score")
  plt %+% pt
  expect_length(plt$tracks[[1]]$elements, 0L)
  expect_length(plt$tracks[[2]]$elements, 1L)
})

test_that("%+% errors when SeqPlot has no tracks", {
  plt <- SeqPlot(windows = .make_windows())
  pt  <- SeqPoint(.make_gr(), yCol = "score")
  expect_error(plt %+% pt)
})

test_that("%+% errors when right-hand side is not a SeqElement", {
  plt <- SeqPlot(windows = .make_windows()) %|% SeqTrack()
  expect_error(plt %+% "not an element")
  expect_error(plt %+% 42)
})


# ── %+% on SeqTrack ──────────────────────────────────────────────────────────

test_that("%+% adds a SeqElement directly to a SeqTrack", {
  trk <- SeqTrack()
  pt  <- SeqPoint(.make_gr(), yCol = "score")
  result <- trk %+% pt
  expect_identical(result, trk)
  expect_length(trk$elements, 1L)
})

test_that("%+% can add multiple elements to a SeqTrack", {
  trk <- SeqTrack()
  gr  <- .make_gr()
  trk %+% SeqPoint(gr, yCol = "score") %+% SeqLine(gr, yCol = "score")
  expect_length(trk$elements, 2L)
})

test_that("%+% errors on SeqTrack when right-hand side is not a SeqElement", {
  trk <- SeqTrack()
  expect_error(trk %+% "bad")
})


# ── Chained operator pipeline ─────────────────────────────────────────────────

test_that("full pipeline SeqPlot() %|% SeqTrack() %+% SeqPoint() builds correctly", {
  gr  <- .make_gr()
  plt <- SeqPlot(windows = .make_windows()) %|%
    SeqTrack() %+%
    SeqPoint(gr, yCol = "score")

  expect_length(plt$tracks, 1L)
  expect_length(plt$tracks[[1]]$elements, 1L)
  expect_true(inherits(plt$tracks[[1]]$elements[[1]], "SeqPoint"))
})

test_that("multi-track pipeline builds the expected structure", {
  gr  <- .make_gr()
  plt <- SeqPlot(windows = .make_windows()) %|%
    SeqTrack() %+% SeqPoint(gr, yCol = "score") %|%
    SeqTrack() %+% SeqLine(gr,  yCol = "score")

  expect_length(plt$tracks, 2L)
  expect_true(inherits(plt$tracks[[1]]$elements[[1]], "SeqPoint"))
  expect_true(inherits(plt$tracks[[2]]$elements[[1]], "SeqLine"))
})
