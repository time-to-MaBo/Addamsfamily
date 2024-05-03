#' row-specific maximum
#'
#' @param x numeric vector.
#' @param ... further arguments passed to max (na.rm)
#'
#' @return logical vector of length dim(x)[1]
#' @export
#'
#' @examples rowall(cbind(c(1,2,3),c(1,1,10))==1)
#' @examples rowall(cbind(c(NA,2,3),c(1,1,10)==1), na.rm = TRUE)
rowall <- function(x, ...){
  x <- as.matrix(x)
  out <- apply(X = x, MARGIN = 1, FUN = all, ...)
  return(out)
}
