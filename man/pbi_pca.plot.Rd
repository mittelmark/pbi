\name{pbi_pca.plot} 
\alias{pbi$pca.plot} 
\alias{pbi_pca.plot} 
\title{Improved bar or screeplot for pca objects.} 
\description{ 
  The function _pca.plot_ provides an improved bar- or screeplot for 
  visualizing the variances of the individual components of an analyses  
  using the function _prcomp_. In contrast to the default plot function  
  this plot visualize cumulative and individual variances in percent. 
} 
\usage{pbi_pca.plot(pca,n=10,type="bar", cex=1.5, legend=TRUE,xlab="Components",ylab="Variance (\%)",pc.col=c("light blue","grey"),...)} 
\arguments{ 
  \item{pca}{ 
    pca object which was created using the function _prcomp_. 
  } 
  \item{n}{ 
    maximal number of components to visualize, default: 10 
  } 
  \item{type}{ 
    plotting type either "bar" or "scree", default: "bar" 
  } 
  \item{cex}{ 
    character expansion for the legend and the screeplot plotting characters, default: 1.5 
  } 
  \item{legend}{ 
    should the legend be displayed on top, default: TRUE 
  } 
  \item{xlab}{ 
     
  } 
  \item{ylab}{ 
     
  } 
  \item{pc.col}{ 
    colors for the PC variances, first individual, second color for the cumulative variance, default: c("light blue","grey") 
  } 
  \item{...}{ 
    additional arguments delegated to the standard plot function 
  } 
} 
\examples{ 
  data(iris) 
  par(mfrow=c(1,2)) 
  pcai=prcomp(iris[,1:4],scale=TRUE) 
  pbi_pca.plot(pcai) 
  pbi_pca.plot(pcai,type="scree",legend=FALSE) 
} 
