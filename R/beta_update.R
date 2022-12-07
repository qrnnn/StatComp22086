#' @title Update the vector.
#' @name beta_update
#' @description Use coordinate descent method to update the vector. The formula is derived from the derivation of single variable problem and the discussion of simple classification. The output is the vector that has been updated.
#' @param x The predictor variable in the elastic net(E-N) problem.
#' @param y The response variable in the E-N problem.
#' @param beta The initial value of the vector to update.
#' @param lambda The coefficient of the penalty term in the E-N problem.
#' @param lambda0 The coefficient of the penalty term L2 in the E-N problem.
#' @export
beta_update <- function(x, y, beta, lambda, lambda0){    
  beta0 <- matrix(1, nrow(beta), 1)
  while(max(abs(beta - beta0)) > 1e-3){
    beta0 <- beta
    for (j in 1:nrow(beta)){   
      temp <- rep(0,nrow(y))
      for (i in 1:nrow(y)){
        temp[i] <- y[i] - sum(beta*x[,i])+beta[j]*x[j,i]
      }
      ru <- sum(x[j,]*temp)
      z <- sum((x^2)[j,]) + 2*lambda0
      if (ru < (-lambda/2)){  
        beta[j] <- ((ru + (lambda/2))/z)
      }else if((ru > (lambda/2))){
        beta[j] <- ((ru - (lambda/2))/z)
      }else{
        beta[j] <- 0
      }    
    }
  }
  return(beta)
}

