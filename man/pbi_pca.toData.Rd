\name{pbi_pca.toData} 
\alias{pbi$pca.toData} 
\alias{pbi_pca.toData} 
\title{Transform prcomp pca object  back to data.} 
\description{ 
  The method allows you transform PCA data back to original data.  
  This can be as well used to eliminate some components and then create 
  data by removing the effect of these components. 
} 
\usage{pbi_pca.toData(pca)} 
\arguments{ 
  \item{pca}{ 
    pca object of class prcomp. 
  } 
} 
\value{return data matrix with the original data.} 
\examples{ 
  pca=prcomp(iris[,1:4]) 
  head(iris[,1:4]) 
  head(pbi_pca.toData(pca)) 
  # remove effect of first component 
  pca$rotation[,1]=0 
  head(pbi_pca.toData(pca)) 
} 
