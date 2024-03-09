\name{pbi_pca.variances} 
\alias{pbi$pca.variances} 
\alias{pbi_pca.variances} 
\title{Return the absolute variance contributions for each variable to each component.} 
\description{ 
  The function eturn the absolute variance contributions for each variable to each component. 
  Every squared loading value for each component and variable  is multiplied 
  with the component importance. The sum of the returned matrix is therefor 1. 
} 
\usage{pbi_pca.variances(pca)} 
\arguments{ 
  \item{pca}{ 
    a PCA object created with `prcomp`. 
  } 
} 
\value{return matrix with absolute variances for each component and variable.} 
\examples{ 
pca=prcomp(USArrests,scale=TRUE) 
round(pbi_pca.variances(pca),2) 
sum(pbi_pca.variances(pca)) 
} 
\seealso{  See also [pbi](#pbi), [pbi$pca.varplot](#pca.varplot) } 
