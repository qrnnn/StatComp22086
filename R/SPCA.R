#' @title SPCA
#' @name SPCA
#' @description PCA is widely used in data processing and dimensionality reduction. However, because the principal component obtained by the method is a linear combination of all the predictive variables, it is difficult to use such principal components to explain the experimental results in many cases. We wanted to filter out the more important (explained with a larger proportion of variance) variables. We need to make the selected principal component as sparse as possible while ensuring that the proportion of variance that can be explained by it is large enough. The SPCA can implement this idea. The idea of SPCA is to transform PCA into a regression-type optimization problem and obtain sparse coefficient matrix by applying lasso constraint to regression coefficient. The output are three items. Loadings: the Loading obtained by SPCA method. Pev: The proportion of variance explained by the principal component. Cumsum_pev: The cumulative proportion of variance explained by the principal component.  
#' @param cormat The correlation coefficient matrix of the original data (if only the original data, then R function cor can be used).
#' @param k The number of principal components to select.
#' @param eps A small value used to control convergence.
#' @param lambda L2 penalty term coefficient.
#' @param lambda1 L1 penalty term coefficient vector.
#' @examples
#' \dontrun{
#' library(elasticnet)
#' data(pitprops)
#' SPCA(pitprops,6,1e-6,0,c(0.06,0.16,0.1,0.5,0.5,0.5))
#' }
#' @import Rcpp elasticnet kableExtra knitr boot bootstrap DAAG MASS microbenchmark purrr
#' @export
SPCA <- function(cormat, k, eps, lambda, lambda1){
  x.eigen <- eigen(cormat)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  x <- v %*% diag(sqrt(d)) %*% t(v)
  n <- nrow(x)
  p <- ncol(x)
  V <- svd(x)$v
  alpha <- V[,1:k]    
  
  Beta_now <- matrix(0, nrow=p, ncol=k)
  Beta_new <- matrix(1.1, nrow=p, ncol=k)
  iter <- 0
  while(max(abs((Beta_new)-(Beta_now)))>eps){
    iter <- iter+1
    Beta_now <- Beta_new
    for (j in 1:k){
      y <- drop(x%*%alpha[,j])
      Beta_new[,j] <- beta_update(as.matrix(x), as.matrix(y), as.matrix(Beta_new[,j]), lambda1[j], lambda)
    }
    tempp <- svd(cormat%*%Beta_new)
    alpha <- (tempp$u)%*%t(tempp$v)    
  }
  normbeta <- sqrt(apply(Beta_new^2, 2, sum))
  normbeta[normbeta==0] <- 1
  beta <- t(t(Beta_new)/normbeta)   
  u <- x%*%beta
  R <- qr.R(qr(u))
  
  pev <- diag(R^2)/sum((svd(x)$d)^2)
  vn <- dimnames(x)[[2]]
  obj <- list(loadings = beta, pev = pev,cumsum_pev = cumsum(pev))
  return(obj)
}
