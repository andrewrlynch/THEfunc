test_that("aes() returns a SeqAes object", {
  a <- aes(y = cell_type)
  expect_s3_class(a, "SeqAes")
})

test_that("aes() deparsed bare symbols to character strings", {
  a <- aes(y = cell_type, fill = mutation_count)
  expect_equal(a[["y"]], "cell_type")
  expect_equal(a[["fill"]], "mutation_count")
})

test_that("aes() with quoted strings are passed through", {
  a <- aes(y = "cell_type")
  # deparse("cell_type") -> '"cell_type"' with outer quotes -- check the actual
  # raw string stored matches the deparsed form of the quoted literal
  expect_equal(a[["y"]], '"cell_type"')
})

test_that("aes() normalises 'col' to 'color' and removes 'col'", {
  a <- aes(col = my_col)
  expect_equal(a[["color"]], "my_col")
  expect_null(a[["col"]])
})

test_that("aes() normalises 'colour' to 'color' and removes 'colour'", {
  a <- aes(colour = my_colour)
  expect_equal(a[["color"]], "my_colour")
  expect_null(a[["colour"]])
})

test_that("aes() 'color' takes precedence over 'col' when both supplied", {
  # Only col/colour → color alias when 'color' is NOT already present
  a <- aes(color = explicit_col, col = other_col)
  # 'color' was explicit, so 'col' alias should not override it
  expect_equal(a[["color"]], "explicit_col")
})

test_that("aes() NULL arguments are excluded from the mapping", {
  a <- aes(y = cell_type, fill = NULL)
  expect_equal(a[["y"]], "cell_type")
  expect_null(a[["fill"]])
})

test_that("aes() with no arguments returns an empty SeqAes", {
  a <- aes()
  expect_s3_class(a, "SeqAes")
  expect_equal(length(a), 0L)
})

test_that("aes() supports all documented aesthetic arguments", {
  a <- aes(x = pos, y = group, fill = score, alpha = conf, size = weight, shape = type)
  expect_equal(a[["x"]], "pos")
  expect_equal(a[["y"]], "group")
  expect_equal(a[["fill"]], "score")
  expect_equal(a[["alpha"]], "conf")
  expect_equal(a[["size"]], "weight")
  expect_equal(a[["shape"]], "type")
})
