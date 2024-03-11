\name{pbi_sem} 
\alias{pbi$sem} 
\alias{pbi_sem} 
\title{standard error of the mean} 
\description{ 
    Calculates the standard error of the mean for a given numerical vector. 
} 
\usage{pbi_sem(x,na.rm=FALSE)} 
\arguments{ 
  \item{x}{ 
    a numerical vector. 
  } 
  \item{na.rm}{ 
    logical vecor indicating if missing values should be removed, default: FALSE 
  } 
} 
\value{return computed standard error of the mean.} 
\examples{ 
  pbi_sem(rnorm(50,mean=10,sd=3)) 
  pbi_sem(rnorm(1000,mean=10,sd=3)) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_cv]{pbi_cv}} } 
