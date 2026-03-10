# ── aes() ─────────────────────────────────────────────────────────────────────

#' Map GRanges metadata columns to visual and positional properties
#'
#' Uses non-standard evaluation so column names can be passed as bare symbols
#' (no quotes needed), matching ggplot2 ergonomics.
#'
#' @param x Bare column name for x-axis position mapping, or `gr` to use the
#'   element's own GRanges positions.
#' @param y Bare column name for y-axis position mapping. Auto-detects scale
#'   type: factor/character columns become discrete, numeric become continuous.
#' @param col,color,colour Bare column name to map to point/line color.
#' @param fill  Bare column name to map to fill color.
#' @param alpha Bare column name to map to transparency (0–1).
#' @param size  Bare column name to map to size.
#' @param shape Bare column name to map to point shape.
#' @return A `SeqAes` object (named list of column-name strings).
#' @export
aes <- function(x = NULL, y = NULL,
                col = NULL, color = NULL, colour = NULL, fill = NULL,
                alpha = NULL, size = NULL, shape = NULL) {
  args    <- as.list(match.call())[-1]
  mapping <- lapply(args, deparse)   # bare symbol CN → character "CN"

  # Normalise col / colour → color
  if ("col"    %in% names(mapping) && !"color" %in% names(mapping))
    mapping[["color"]] <- mapping[["col"]]
  if ("colour" %in% names(mapping) && !"color" %in% names(mapping))
    mapping[["color"]] <- mapping[["colour"]]
  mapping[c("col", "colour")] <- NULL

  structure(mapping, class = "SeqAes")
}

# ── Internal aesthetic resolver ────────────────────────────────────────────────

#' Resolve a SeqAes + SeqScale into concrete per-observation vectors
#'
#' @param data_mcols  data.frame of GRanges metadata columns.
#' @param aes_obj     A `SeqAes` object (from `aes()`), or `NULL`.
#' @param scale_obj   A `SeqScale` object (from `seq_scale_*()`), or `NULL`.
#' @param n           Number of observations.
#' @param default_color Fallback color when no mapping is specified.
#' @return Named list with resolved vectors for `color`, `fill`, `alpha`,
#'   `size`, `shape` (only those present in `aes_obj`).
#'
#' @keywords internal
.resolve_aes <- function(data_mcols, aes_obj, scale_obj, n,
                          default_color = "#1C1B1A") {
  result <- list(color = rep(default_color, n))
  if (is.null(aes_obj)) return(result)

  for (aes_name in c("color", "fill", "alpha", "size", "shape")) {
    col_name <- aes_obj[[aes_name]]
    if (is.null(col_name) || !col_name %in% names(data_mcols)) next
    raw <- data_mcols[[col_name]]
    sc  <- if (!is.null(scale_obj) && scale_obj$aesthetic == aes_name) scale_obj else NULL

    if (aes_name %in% c("color", "fill")) {
      is_discrete <- is.factor(raw) || is.character(raw) ||
                     (!is.null(sc) && sc$type == "discrete")

      if (is_discrete) {
        raw_f  <- as.factor(raw)
        lvls   <- levels(raw_f)
        pal    <- if (!is.null(sc) && !is.null(sc$values)) {
          sc$values
        } else if (!is.null(sc) && !is.null(sc$palette)) {
          sc$palette(length(lvls))
        } else {
          flexoki_palette(length(lvls))
        }
        names(pal)              <- lvls[seq_along(pal)]
        resolved                <- pal[as.character(raw_f)]
        resolved[is.na(resolved)] <- if (!is.null(sc)) sc$na_value %||% "grey80" else "grey80"

      } else {
        raw_n  <- as.numeric(raw)
        lims   <- if (!is.null(sc) && !is.null(sc$limits)) sc$limits
                  else range(raw_n, na.rm = TRUE)
        pal_nm <- if (!is.null(sc)) sc$palette %||% "viridis" else "viridis"
        stops  <- switch(pal_nm,
                    viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
                    plasma  = c("#0d0887", "#cc4778", "#f0f921"),
                    magma   = c("#000004", "#b63679", "#fcfdbf"),
                    blues   = c("#f7fbff", "#2171b5", "#08306b"),
                    reds    = c("#fff5f0", "#ef3b2c", "#67000d"),
                    stop("Unknown palette: ", pal_nm))
        ramp   <- grDevices::colorRamp(stops)
        t_vals <- pmax(0, pmin(1, (raw_n - lims[1]) / diff(lims)))
        cols   <- ramp(t_vals)
        resolved <- grDevices::rgb(cols[, 1], cols[, 2], cols[, 3],
                                   maxColorValue = 255)
        resolved[is.na(raw_n)] <- if (!is.null(sc)) sc$na_value %||% "grey80" else "grey80"
      }

      result[[aes_name]] <- unname(resolved)

    } else {
      # size, alpha, shape — numeric pass-through
      result[[aes_name]] <- as.numeric(raw)
    }
  }
  result
}

# ── Scale functions ────────────────────────────────────────────────────────────

#' Continuous color scale for SeqAes mappings
#'
#' @param palette One of `"viridis"`, `"plasma"`, `"magma"`, `"blues"`, `"reds"`.
#' @param limits  Optional numeric vector of length 2 clamping the scale range.
#' @param na_value Color for NA values (default `"grey80"`).
#' @param midpoint Optional midpoint for diverging scales (not yet implemented).
#' @return A `SeqScaleContinuous` / `SeqScale` object.
#' @export
seq_scale_color_continuous <- function(palette  = "viridis",
                                       limits   = NULL,
                                       na_value = "grey80",
                                       midpoint = NULL) {
  structure(
    list(aesthetic = "color", type = "continuous",
         palette = palette, limits = limits,
         na_value = na_value, midpoint = midpoint),
    class = c("SeqScaleContinuous", "SeqScale")
  )
}

#' Discrete color scale for SeqAes mappings
#'
#' @param values  Optional named character vector mapping level names to colors.
#' @param palette Optional function `function(n)` returning `n` hex colors.
#'   Falls back to `flexoki_palette()` if `NULL`.
#' @param na_value Color for NA / unmatched values (default `"grey80"`).
#' @return A `SeqScaleDiscrete` / `SeqScale` object.
#' @export
seq_scale_color_discrete <- function(values   = NULL,
                                     palette  = NULL,
                                     na_value = "grey80") {
  structure(
    list(aesthetic = "color", type = "discrete",
         values = values, palette = palette, na_value = na_value),
    class = c("SeqScaleDiscrete", "SeqScale")
  )
}

#' @rdname seq_scale_color_continuous
#' @export
seq_scale_fill_continuous <- function(...) {
  s <- seq_scale_color_continuous(...); s$aesthetic <- "fill"; s
}

#' @rdname seq_scale_color_discrete
#' @export
seq_scale_fill_discrete <- function(...) {
  s <- seq_scale_color_discrete(...); s$aesthetic <- "fill"; s
}



# ── Position scale constructors ──────────────────────────────────────────────

#' Genomic position scale
#'
#' Creates a position scale based on genomic coordinates from a `GRanges` object.
#' Used for axes that represent physical genome positions (bp, kb, Mb).
#'
#' @param windows A `GRanges` object defining the genomic windows.
#' @param scale_factor Numeric vector of per-window scale factors controlling
#'   the unit label (1e-6 = Mb, 1e-3 = kb, 1 = bp). If `NULL`, reads from
#'   `mcols(windows)$scale` or defaults to `1e-6`.
#' @return A `SeqScaleGenomic` / `SeqPositionScale` object.
#' @export
seq_scale_genomic <- function(windows, scale_factor = NULL) {
  stopifnot(inherits(windows, "GRanges"))
  if (is.null(scale_factor)) {
    scale_factor <- if ("scale" %in% names(S4Vectors::mcols(windows)))
      S4Vectors::mcols(windows)$scale
    else
      rep(1e-6, length(windows))
  }
  structure(
    list(type = "genomic", windows = windows, scale_factor = scale_factor),
    class = c("SeqScaleGenomic", "SeqPositionScale")
  )
}

#' Continuous position scale
#'
#' Creates a numeric position scale for axes displaying continuous data.
#'
#' @param limits Optional numeric vector of length 2 clamping the axis range.
#'   If `NULL`, the range is auto-computed from element data.
#' @param n_breaks Target number of pretty breaks (default 5).
#' @return A `SeqScaleContinuous_Pos` / `SeqPositionScale` object.
#' @export
seq_scale_continuous <- function(limits = NULL, n_breaks = 5) {
  structure(
    list(type = "continuous", limits = limits, n_breaks = n_breaks),
    class = c("SeqScaleContinuous_Pos", "SeqPositionScale")
  )
}

#' Discrete position scale
#'
#' Creates a categorical position scale for axes with discrete levels
#' (e.g., cell types, sample names).
#'
#' @param levels Character vector of category levels, in display order.
#'   If `NULL`, levels are auto-detected from element data.
#' @param labels Optional display labels (same length as `levels`).
#'   If `NULL`, level names are used as labels.
#' @return A `SeqScaleDiscrete_Pos` / `SeqPositionScale` object.
#' @export
seq_scale_discrete <- function(levels = NULL, labels = NULL) {
  structure(
    list(type = "discrete", levels = levels, labels = labels),
    class = c("SeqScaleDiscrete_Pos", "SeqPositionScale")
  )
}