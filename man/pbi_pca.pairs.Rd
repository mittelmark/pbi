\name{pbi_pca.pairs} 
\alias{pbi$pca.pairs} 
\alias{pbi_pca.pairs} 
\title{Improved pairs plot for pca objects.} 
\description{ 
  The function _pca.pairs_ provides an improved pairs plot for 
  visualizing the pairwise scores of the individual components of an analyses  
  using the function _prcomp_. In contrast to the default pairs function  
  this plot visualizes in the diagonal as well the variances and  
  a density line for the component scores. 
} 
\usage{pbi_pca.pairs(pca,n=10,groups=NULL,col='black',pch=19,legend=FALSE,...)} 
\arguments{ 
  \item{pca}{ 
    pca object which was created using the function _prcomp_. 
  } 
  \item{n}{ 
    maximal number of components to visualize, default: 10 
  } 
  \item{groups}{ 
    vector with classes having the same length than the input matrix for prcomp has rows, default: NULL 
  } 
  \item{col}{ 
    colors for the plotting, character, default: 'black' 
  } 
  \item{pch}{ 
    plotting, symbol, default: 19 
  } 
  \item{legend}{ 
    should the legend be displayed on top, default: FALSE 
  } 
  \item{...}{ 
    additional arguments delegated to the standard _pairs_ function 
  } 
} 
\examples{ 
  data(iris) 
  pci=prcomp(iris[,1:4],scale=TRUE) 
  pbi_pca.pairs(pci,pch=15,groups=iris[,5], 
     legend=TRUE,oma=c(5,4,4,4),col=as.numeric(iris[,5])+1) 
} 
