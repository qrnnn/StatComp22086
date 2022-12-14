# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title EM of certain case
#' @description EM method of certain case: we only know that Xi's have distribution Exp(a)(i.i.d) and Xi falls in (ui,vi) for i=1,...,n. This function can use to get the unknown parameter a. 
#' @param lambda1 This method requires two initialization parameters, it's the first one.
#' @param lambda2 This method requires two initialization parameters, it's the second one.
#' @param u A vector consisting of the left bounds of several known intervals.
#' @param v A vector consisting of the right bounds of several known intervals.
#' @return The estimate of a from the EM method.
#' @examples
#' \dontrun{
#' lambda1 <- 2
#' lambda2 <- 1
#' u <- c(11,8,27,13,16,0,23,10,24,2)
#' v <- c(12,9,28,14,17,1,24,11,25,3)
#' emc <- EMC(lambda1,lambda2,u,v)
#' emc
#' }
#' @export
EMC <- function(lambda1, lambda2, u, v) {
    .Call('_StatComp22086_EMC', PACKAGE = 'StatComp22086', lambda1, lambda2, u, v)
}

gibbsC <- function(mu1, mu2, sigma1, sigma2, rho, N, burn) {
    .Call('_StatComp22086_gibbsC', PACKAGE = 'StatComp22086', mu1, mu2, sigma1, sigma2, rho, N, burn)
}

