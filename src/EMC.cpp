#include <Rcpp.h>
using namespace Rcpp;

//' @title EM of certain case
//' @description EM method of certain case: we only know that Xi's have distribution Exp(a)(i.i.d) and Xi falls in (ui,vi) for i=1,...,n. This function can use to get the unknown parameter a. 
//' @param lambda1 This method requires two initialization parameters, it's the first one.
//' @param lambda2 This method requires two initialization parameters, it's the second one.
//' @param u A vector consisting of the left bounds of several known intervals.
//' @param v A vector consisting of the right bounds of several known intervals.
//' @return The estimate of a from the EM method.
//' @examples
//' \dontrun{
//' lambda1 <- 2
//' lambda2 <- 1
//' u <- c(11,8,27,13,16,0,23,10,24,2)
//' v <- c(12,9,28,14,17,1,24,11,25,3)
//' emc <- EMC(lambda1,lambda2,u,v)
//' emc
//' }
//' @export
// [[Rcpp::export]]
double EMC(double lambda1, double lambda2, NumericVector u, NumericVector v){
  double eps = pow(10,-5);
  int n = sizeof(u);
  while(abs(lambda1-lambda2)>=eps){
    lambda1 = lambda2;
    NumericVector x = exp(-lambda1*u);
   NumericVector y = exp(-lambda1*v);
    lambda2 = n/sum(((lambda1*u+1)*x-(lambda1*v+1)*y)/(lambda1*(x-y)));
  }
  return(lambda2);
}
