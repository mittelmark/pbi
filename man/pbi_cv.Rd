\name{pbi_cv} 
\alias{pbi$cv} 
\alias{pbi_cv} 
\title{Coefficient of variation} 
\description{ 
  The function \code{pbi_cv} calculates the coefficient of variation, which is 
  an unit-less measure of data scatter. It is calculated as: \emph{(100*sd(x))/mean(x)}. 
  All data values has to be non-negative. 
} 
\usage{pbi_cv(x,na.rm=FALSE)} 
\arguments{ 
  \item{x}{ 
    vector with positive numerical values. 
  } 
  \item{na.rm}{ 
    should NA's be removed, default: FALSE 
  } 
} 
\value{return numerical value for the coefficient of variation.} 
\examples{ 
  pbi_cv(rnorm(20,mean=100,sd=4)) 
  pbi_cv(c(1,2,3,4)) 
  pbi_cv(c(1,2,3,4,NA)) 
  pbi_cv(c(1,2,3,4,NA),na.rm=TRUE) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_sem]{pbi_sem}} } 
