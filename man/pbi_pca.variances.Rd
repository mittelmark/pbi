\name{pbi_pca.variances} 
\alias{pbi$pca.variances} 
\alias{pbi_pca.variances} 
\title{Absolute variance contributions of variables to PC's} 
\description{ 
  The function returns the absolute variance contributions for each variable to each component. 
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
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_pca.biplot]{pbi_pca.biplot}},  
\code{\link[pbi:pbi_pca.corplot]{pbi_pca.corplot}},  
\code{\link[pbi:pbi_pca.pairs]{pbi_pca.pairs}}, 
\code{\link[pbi:pbi_pca.plot]{pbi_pca.plot}}, 
\code{\link[pbi:pbi_pca.toData]{pbi_pca.toData}}, 
\code{\link[pbi:pbi_pca.varplot]{pbi_pca.varplot}} 
} 
