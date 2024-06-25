#' The wrapper of the ns function
#'
#' @param x a predictor matrix
#' @param y a response vector
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h the number of observations to minimise sum of squares over; by default a sequence from
#' 75 to 100 percent of sample size (in increments of 5 percent)
#' @param max.ns.iter the maximum number of neighbourhood search iterations allowed
#' @param max.gd.iter the maximum number of gradient descent iterations allowed per value of
#' \code{k} and \code{h}
#' @param eps a numerical tolerance parameter used to declare convergence
#'
#' @return a list with the following components:
#' \item{beta}{an array of estimated regression coefficients; columns correspond to \code{k} and
#' matrices to \code{h}}
#' \item{eta}{an array of estimated residuals; columns correspond to \code{k} and matrices to
#' \code{h}}
#'
#' @export

heuristics_ns <- function(x, y, k, h, max_ns_iter, max_gd_iter, eps) {
    .Call(`_robustsubsets_ns`, x, y, k, h, max_ns_iter, max_gd_iter, eps)
}
