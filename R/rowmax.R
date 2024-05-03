#' row-specific maximum
#'
#' @param x numeric vector.
#' @param ... further arguments passed to max (na.rm)
#'
#' @return numeric vector of length dim(x)[1] with row-secific maxima
#' @export
#'
#' @examples rowmax(cbind(c(1,2,3),c(1,1,10)))
#' @examples rowmax(cbind(c(NA,2,3),c(1,1,10)), na.rm = TRUE)
rowmax <- function(x, ...){
  x <- as.matrix(x)
  out <- apply(X = x, MARGIN = 1, FUN = max, ...)
  return(out)
}
