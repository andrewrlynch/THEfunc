#' drawSeqArch
#'
#' General function for drawing SeqLink arches
#'
#' @param x0,x1,y0,y1,top0,top1 Canvas coordinates. y0 and y1 are bottoms of arch stems. top0 and top1 are the tops of arch stems
#' @param orientation Orientation of the link. In general, "-" arches will point down, "+" will point up.
#' @param curve Argument to set curve type or degree. Can be one of  "length", "equal", or number to control height of arc (default is 0.2).
#' @param stemWidth,arcWidth,arcColor,stemColor Aesthetic control over arches.
#'
#' @export
drawSeqArch <- function(x0, y0, x1, y1, top0, top1,
                        orientation = "*", curve = "length",
                        stemWidth = 1, arcWidth = 1,
                        arcColor = "black", stemColor = "black") {

  #upward <- orientation %in% c("*", "up", "+", "-/+", "-/-")
  #sign <- ifelse(upward, 1, -1)

  span <- abs(x1 - x0)

  # Determine vertical arc "bulge" height
  if (is.numeric(curve)) {
    curve_offset <- curve
  } else if (curve == "equal") {
    curve_offset <- 0.2  # default uniform curve
  } else if (curve == "length") {
    curve_offset <- span * 0.2  # scale with span (same multiplier)
  } else {
    warning("Unknown curve value in drawSeqArch(); using default 0.2")
    curve_offset <- span * 0.2
  }

  curve_offset1 <- ifelse(grepl("^\\*|^\\+", orientation), curve_offset, curve_offset * -1)
  curve_offset2 <- ifelse(grepl("\\*$|\\+$", orientation), curve_offset, curve_offset * -1)


  # Controls how broad the arch appears (smaller = broader)
  ctrl_spread <- 0

  # Horizontal offset from endpoints
  dx <- abs(x1 - x0)
  ctrl_dx <- dx * ctrl_spread

  mid_y <- max(top0, top1)
  ctrl_x1 <- x0 + ctrl_dx
  ctrl_x2 <- x1 - ctrl_dx

  ctrl_y1 <- mid_y + curve_offset1
  ctrl_y2 <- mid_y + curve_offset2

  P0 <- c(x0, top0)
  P1 <- c(ctrl_x1, ctrl_y1)
  P2 <- c(ctrl_x2, ctrl_y2)
  P3 <- c(x1, top1)

  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^3 * P0[1] + 3 * (1 - t)^2 * t * P1[1] +
    3 * (1 - t) * t^2 * P2[1] + t^3 * P3[1]
  bez_y <- (1 - t)^3 * P0[2] + 3 * (1 - t)^2 * t * P1[2] +
    3 * (1 - t) * t^2 * P2[2] + t^3 * P3[2]

  grid.lines(bez_x, bez_y, gp = gpar(col = arcColor, lwd = arcWidth))

  # Draw vertical stems
  grid.segments(x0, y0, x0, top0, gp = gpar(col = stemColor, lwd = stemWidth))
  grid.segments(x1, y1, x1, top1, gp = gpar(col = stemColor, lwd = stemWidth))
}


drawSeqString <- function(x0, y0, x1, y1,
                          strand1 = "*",
                          strand2 = "*",
                          orientation = "*",
                          type = c("c", "s"),
                          bulge = 0.04,
                          lwd = 1.5,
                          col = "red",
                          alpha = 1) {

  if (length(x0)==0 || length(y0)==0 || length(x1)==0 || length(y1)==0)
    return(invisible(NULL))

  x0 <- as.numeric(x0)[1]; y0 <- as.numeric(y0)[1]
  x1 <- as.numeric(x1)[1]; y1 <- as.numeric(y1)[1]
  if (!is.finite(x0) || !is.finite(y0) || !is.finite(x1) || !is.finite(y1))
    return(invisible(NULL))

  strand1 <- as.character(strand1[1])
  strand2 <- as.character(strand2[1])
  orientation <- tolower(as.character(orientation[1]))

  type <- tolower(as.character(type[1]))
  if (!type %in% c("c","s")) type <- "c"

  # keep stable left->right; if we swap endpoints, swap strand roles too
  if (x1 < x0) {
    tx <- x0; x0 <- x1; x1 <- tx
    ty <- y0; y0 <- y1; y1 <- ty
    ts <- strand1; strand1 <- strand2; strand2 <- ts
  }

  dx <- x1 - x0
  dy <- y1 - y0
  if (!is.finite(dx) || dx <= 0) return(invisible(NULL))

  bx <- as.numeric(bulge)[1]
  if (!is.finite(bx)) bx <- 0.07
  bx <- max(0, min(0.35, bx))   # npc plot-scale clamp

  # direction (right=+1, left=-1)
  dir <- 1

  # explicit orientation wins if provided
  if (orientation %in% c("-", "down", "left", "neg")) {
    dir <- -1
  } else if (orientation %in% c("+", "up", "right", "pos")) {
    dir <- 1
  } else {
    # otherwise derive from strand1 when available
    if (strand1 == "-") dir <- -1
    if (strand1 == "+") dir <- 1
  }

  # y-control points
  #yP1 <- y0 + dy/3
  #yP2 <- y0 + 2*dy/3
  yP1 <- y0
  yP2 <- y1

  # x-control points: horizontal bulge
  same_strand <- (strand1 %in% c("+","-")) && (strand2 %in% c("+","-")) && (strand1 == strand2)

  if (same_strand) {
    # C curve: both controls bend same direction
    xP1 <- x0 + dir * bx
    xP2 <- x1 + dir * bx
  } else {
    # S curve: controls bend opposite directions
    xP1 <- x0 + dir * bx
    xP2 <- x1 - dir * bx
  }

  xP1 <- max(0, min(1, xP1))
  xP2 <- max(0, min(1, xP2))

  gp <- grid::gpar(col = grDevices::adjustcolor(col, alpha.f = alpha), lwd = lwd)

  grid::grid.draw(grid::bezierGrob(
    x = c(x0, xP1, xP2, x1),
    y = c(y0, yP1, yP2, y1),
    default.units = "npc",
    gp = gp
  ))

  invisible(NULL)
}