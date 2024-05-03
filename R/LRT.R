#' LRT
#'
#' @param l_h1 loglihood value of larger model
#' @param l_h0 loglihood value of nested model
#' @param df degrees of fredom
#'
#' @return numeric vector with LRT, df, -value
#' @export
#'
#' @examples
LRT <- function(l_h1, l_h0, df){

  lrt <- 2*(l_h1 - l_h0)
  p_value <- stats::pchisq(q = lrt, df = df, lower.tail = FALSE)

  OUT <- matrix(c(lrt, df, p_value), ncol = 3)
  colnames(OUT) <- c('LRT', 'df', 'p_value')

  return(OUT)
}
