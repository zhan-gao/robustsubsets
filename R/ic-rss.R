#' Information Criteria to choose k and h
#'
#' @param rss_obj the return object from the function rss

#' @return The chosen k and h
#'
#' @export
#'


ic_rss <- function(rss_obj) {
   n <- dim(rss_obj$eta)[1]
   h <- rss_obj$h
   nh <- length(h)
   k <- rss_obj$k
   nk <- length(k)

   # degree of freedom
   df_mat <- (k %*% t(rep(1, nh))) + (rep(1, nk) %*% t(as.matrix(n - h))) + 1 # intercept
   # information criteria
   ic_mat <- n * log(rss_obj$objval / n) + df_mat * log(n)
   pos <- which(ic_mat == min(ic_mat), arr.ind = TRUE)
   k_out <- k[pos[1]]
   h_out <- h[pos[2]]
   return(
       list(k = k_out,
            h = h_out)
   )
}



