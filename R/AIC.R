#' AIC
#'
#' @param loglihood loglihood value
#' @param num_par number of parameters
#'
#' @return numeric. AIC. sclar
#' @export
#'
#' @examples
aic <- function(loglihood, num_par){

  aic <- 2*(num_par - loglihood)
  return(aic)
}
