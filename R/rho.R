#' Convex function for estimation of M-estimator
#'
#' @param u Value of variable.
#' @param func If 1 (default), the convex function is absolute function. If 2, the convex function is square function.
#'
#' @return value of convex function
#'
#' @export
rho <- function(u, func=1){
if(func==1){
  return(abs(u))
}else if(func==2){
  return(u^2)
}else{
  return(u)
}
}
