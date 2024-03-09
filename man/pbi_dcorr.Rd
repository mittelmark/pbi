\name{pbi_dcorr} 
\alias{pbi$dcorr} 
\alias{pbi_dcorr} 
\title{Calculate pairwise correlations for a given data frame or matrix including their p-values.} 
\description{ 
  The function is an extension to the standard _stats::cor_ function, it calculates as well 
  the p-values for the pairwise aoosciations and returns them in a matrix as well. 
} 
\usage{pbi_dcorr(data,method='pearson',use='pairwise.complete.ob')} 
\arguments{ 
  \item{data}{ 
    matrix or data frame where the variables are in the columns, NA's are allowed. 
  } 
  \item{method}{ 
    type of correlation to be determined, either 'pearson', 'spearman' or 'kendall', default: 'pearson' 
  } 
  \item{use}{ 
    how to deal with NA's, default: 'pairwise.complete.obs' 
  } 
} 
\value{return list with the following components: 
 
> - _estimate_ - matrix with correlation values 
  - _p.value_ - matrix with p-values 
  - _method_ - character string with the used correlation method 
  } 
\examples{ 
  data(swiss) 
  lapply(pbi_dcorr(swiss)[1:2],round,2) 
} 
\seealso{  See also [pbi](#home), [pbi_dcorrplot](#dcorrplot) } 
