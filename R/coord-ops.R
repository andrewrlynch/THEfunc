.SeqUnitWindow <- function() {
  # Use a valid IRanges (>=1). We won't use genomic scaling for legend anyway.
  win <- GenomicRanges::GRanges("legend", IRanges::IRanges(1, 2))
  S4Vectors::mcols(win)$scale <- 1  # prevent Mb scaling surprises (even though axes are off)
  win
}

#' clipToXscale
#'
#' Clipping for ranges intersected by x axis limits
#'
#' @param x0,x1,xscale Coordinates for clipping
#'
#' @export
clipToXscale <- function(x0, x1, xscale) {
  keep <- x1 >= xscale[1] & x0 <= xscale[2]
  x0_clipped <- pmax(x0[keep], xscale[1])
  x1_clipped <- pmin(x1[keep], xscale[2])
  return(list(x0 = x0_clipped, x1 = x1_clipped, mask = keep))
}



#' convertDataToGrid
#'
#' Converts data coordinates to grid coordinates
#'
#' @param x,y,xscale,yscale Coordinates for clipping
#'
#' @export
convertDataToGrid <- function(x, y, xscale, yscale) {
  u <- (x - xscale[1]) / diff(xscale)
  v <- (y - yscale[1]) / diff(yscale)
  c(u, v)
}



#' gridToCanvas
#'
#' Converts grid coordinates to canvas coordinates
#'
#' @param u,v,panel_meta Coordinates for clipping
#'
#' @export
gridToCanvas <- function(u, v, panel_meta) {
  x <- panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0)
  y <- panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)
  c(x, y)
}