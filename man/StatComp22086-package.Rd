\name{StatComp22086-package}
\alias{StatComp22086-package}
\alias{StatComp22086}
\docType{package}
\title{
  EM for a certain case and SPCA
}
\description{
  In addition to homework, this package includes three functions. Two of them make up the SPCA method. This method can make the loadings of the selected principal component as sparse as possible under the condition that the proportion of the interpretable variance of the selected principal component is large enough. Another function is the EM algorithm written in RCPP for a particular case.
}
\details{
  beta_update: Update the parameters (Coordinate descent method).
  SPCA: The function of SPCA method. Output parameters: loadings(loading obtained by SPCA method), pev(the proportion of variance explained by the principal component), cumsum_pev(the cumulative proportion of variance explained by the principal component).
  EMC: EM method of certain case. We only know that Xi's have distribution Exp(a)(i.i.d) and Xi falls in (ui,vi) for i=1,...,n. This function can use to get the unknown parameter a. 
}
\author{
  Qi Ruonan
  Maintainer: Qi Ruonan <qrn@mail.ustc.edu.cn>
}
\keyword{ package }
\examples{
  \dontrun{
    library(elasticnet)
    data(pitprops)
    SPCA(pitprops,6,1e-6,0,c(0.06,0.16,0.1,0.5,0.5,0.5))
  }
}
