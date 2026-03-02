# Tests for the internal .resolve_aes() function.
# Access via THEfunc:::.resolve_aes() since it is not exported.

.resolve <- THEfunc:::.resolve_aes


# ── NULL / no-op cases ────────────────────────────────────────────────────────

test_that(".resolve_aes() with NULL aes returns default color vector", {
  df     <- data.frame(x = 1:4)
  result <- .resolve(df, NULL, NULL, n = 4, default_color = "#ABCDEF")
  expect_equal(result$color, rep("#ABCDEF", 4))
})

test_that(".resolve_aes() with NULL aes returns exactly n entries", {
  df     <- data.frame(x = 1:7)
  result <- .resolve(df, NULL, NULL, n = 7, default_color = "grey")
  expect_length(result$color, 7L)
})


# ── Discrete color mapping ────────────────────────────────────────────────────

test_that(".resolve_aes() maps character color column to per-obs colors", {
  df <- data.frame(group = c("A", "B", "A"), stringsAsFactors = FALSE)
  a  <- aes(color = group)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  expect_length(result$color, 3L)
  # Observations with the same level get the same color
  expect_equal(result$color[1], result$color[3])
  # Different levels get different colors
  expect_false(result$color[1] == result$color[2])
})

test_that(".resolve_aes() maps factor color column correctly", {
  df <- data.frame(group = factor(c("X", "Y", "X"), levels = c("X", "Y")))
  a  <- aes(color = group)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  expect_equal(result$color[1], result$color[3])  # both "X"
  expect_false(result$color[1] == result$color[2])
})

test_that(".resolve_aes() uses named values from SeqScaleDiscrete", {
  df  <- data.frame(cat = c("A", "B", "A"), stringsAsFactors = FALSE)
  a   <- aes(color = cat)
  sc  <- seq_scale_color_discrete(values = c(A = "#FF0000", B = "#0000FF"))
  result <- .resolve(df, a, sc, n = 3, default_color = "grey")
  expect_equal(result$color[1], "#FF0000")  # A
  expect_equal(result$color[2], "#0000FF")  # B
  expect_equal(result$color[3], "#FF0000")  # A again
})

test_that(".resolve_aes() uses palette function from SeqScaleDiscrete", {
  df   <- data.frame(cat = c("P", "Q"), stringsAsFactors = FALSE)
  a    <- aes(color = cat)
  pal  <- function(n) c("#111111", "#222222")[seq_len(n)]
  sc   <- seq_scale_color_discrete(palette = pal)
  result <- .resolve(df, a, sc, n = 2, default_color = "grey")
  expect_equal(result$color[1], "#111111")
  expect_equal(result$color[2], "#222222")
})

test_that(".resolve_aes() applies na_value for unmatched discrete levels", {
  df  <- data.frame(cat = c("A", NA_character_), stringsAsFactors = FALSE)
  a   <- aes(color = cat)
  sc  <- seq_scale_color_discrete(values = c(A = "#FF0000"), na_value = "#CCCCCC")
  result <- .resolve(df, a, sc, n = 2, default_color = "grey")
  expect_equal(result$color[2], "#CCCCCC")
})


# ── Continuous color mapping ──────────────────────────────────────────────────

test_that(".resolve_aes() maps numeric column to hex colors (viridis default)", {
  df <- data.frame(val = c(0, 50, 100))
  a  <- aes(color = val)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  # All results should be valid hex colors
  expect_match(result$color[1], "^#[0-9A-Fa-f]{6}$")
  expect_match(result$color[2], "^#[0-9A-Fa-f]{6}$")
  expect_match(result$color[3], "^#[0-9A-Fa-f]{6}$")
  # Extremes should differ
  expect_false(result$color[1] == result$color[3])
})

test_that(".resolve_aes() viridis min value matches first stop", {
  df <- data.frame(val = c(0, 100))
  a  <- aes(color = val)
  sc <- seq_scale_color_continuous(palette = "viridis", limits = c(0, 100))
  result <- .resolve(df, a, sc, n = 2, default_color = "grey")
  expect_equal(toupper(result$color[1]), "#440154")
})

test_that(".resolve_aes() viridis max value matches last stop", {
  df <- data.frame(val = c(0, 100))
  a  <- aes(color = val)
  sc <- seq_scale_color_continuous(palette = "viridis", limits = c(0, 100))
  result <- .resolve(df, a, sc, n = 2, default_color = "grey")
  expect_equal(toupper(result$color[2]), "#FDE725")
})

test_that(".resolve_aes() uses custom limits to clamp mapping", {
  df <- data.frame(val = c(-10, 0, 200))
  a  <- aes(color = val)
  sc <- seq_scale_color_continuous(palette = "viridis", limits = c(0, 100))
  result <- .resolve(df, a, sc, n = 3, default_color = "grey")
  # -10 clamped to 0 → same color as 0
  expect_equal(result$color[1], result$color[2])
  # 200 clamped to 100 → same color as the full-range max
  sc2    <- seq_scale_color_continuous(palette = "viridis", limits = c(0, 100))
  ref    <- .resolve(data.frame(val = 100), aes(color = val), sc2, n = 1, "grey")
  expect_equal(result$color[3], ref$color[1])
})

test_that(".resolve_aes() applies na_value for NA continuous entries", {
  df <- data.frame(val = c(0, NA, 100))
  a  <- aes(color = val)
  sc <- seq_scale_color_continuous(na_value = "#AABBCC")
  result <- .resolve(df, a, sc, n = 3, default_color = "grey")
  expect_equal(result$color[2], "#AABBCC")
})


# ── Fill mapping ──────────────────────────────────────────────────────────────

test_that(".resolve_aes() resolves fill separately from color", {
  df <- data.frame(grp = c("X", "Y"), val = c(0.0, 1.0),
                   stringsAsFactors = FALSE)
  a  <- aes(color = grp, fill = val)
  sc <- seq_scale_fill_continuous(palette = "plasma", limits = c(0, 1))
  result <- .resolve(df, a, sc, n = 2, default_color = "grey")
  # fill should be resolved
  expect_match(result$fill[1], "^#[0-9A-Fa-f]{6}$")
  # color was discrete and also resolved
  expect_false(result$color[1] == result$color[2])
})


# ── Passthrough aesthetics (size, alpha, shape) ───────────────────────────────

test_that(".resolve_aes() passes numeric column through for 'size'", {
  df <- data.frame(sz = c(1.0, 2.5, 0.5))
  a  <- aes(size = sz)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  expect_equal(result$size, c(1.0, 2.5, 0.5))
})

test_that(".resolve_aes() passes numeric column through for 'alpha'", {
  df <- data.frame(alpha_col = c(0.2, 0.8, 1.0))
  a  <- aes(alpha = alpha_col)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  expect_equal(result$alpha, c(0.2, 0.8, 1.0))
})

test_that(".resolve_aes() only resolves aesthetics present in aes()", {
  df     <- data.frame(val = 1:3)
  a      <- aes(color = val)
  result <- .resolve(df, a, NULL, n = 3, default_color = "grey")
  # 'fill', 'alpha', 'size', 'shape' should not be added
  expect_null(result$fill)
  expect_null(result$size)
  expect_null(result$alpha)
  expect_null(result$shape)
})
