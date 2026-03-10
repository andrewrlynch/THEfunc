#' Load colors from Kepano's Flexoki palette
#'
#' Loads colors from Kepano's Flexoki palette.
#'
#' @param n Number of colors to be produced. If > 9, an equally distributed color ramp will be produced (and probably look bad).
#'
#' @return Vector of hex colors
#' @export
flexoki_palette <- function(n, shade = 400) {
  flexoki600 <- c("#AF3029","#205EA6","#66800B","#BC5215","#24837B","#5E409D","#AD8301","#A02C6D","#6F6E69")
  flexoki400 <- c("#D14D41","#4385BE","#879A39","#DA702C","#3AA99F","#8B7EC8","#D0A215","#CE5D97","#9F9D96")
  flexoki300 <- c("#E8705F","#66A0C8","#A0AF54","#EC8B49","#5ABDAC","#A699D0","#DFB431","#E47DA8","#B7B5AC")
  flexoki200 <- c("#F89A8A","#92BFDB","#BEC97E","#F9AE77","#87D3C3","#C4B9E0","#ECCB60","#F4A4C2","#CECDC3")

  base <- if(shade == 600){flexoki600} else if (shade == 400){flexoki400} else if (shade == 300){flexoki300} else if (shade == 200){flexoki200} else {stop("shade must be 200, 300, 400, or 600")}

  if (n <= length(base)) {
    return(base[seq_len(n)])
  } else {
    # Interpolation in HCL space
    colorRampPalette(base, space = "Lab")(n)
  }
}