\name{pbi_pca.varplot} 
\alias{pbi$pca.varplot} 
\alias{pbi_pca.varplot} 
\title{PCA variance plot to total variances for each component and variable.} 
\description{ 
    The function provides a PCA matrix plot to show associations  
  between PCs and variables. Sown are the squared values, but retaining the  
  original sign of the the variances. So the abolute sum of all values should be one. 
} 
\usage{pbi_pca.varplot(pca,pcs=10,main="Variance plot",cex.lab=1.5,cex.sym=8, cex.var=1, pch=16,...)} 
\arguments{ 
  \item{pca}{ 
    pca object which was created using the function _prcomp_. 
  } 
  \item{pcs}{ 
    number of the first PC's or vector of PC names to be plotted, default: 10 (or max number of PC's) 
  } 
  \item{main}{ 
    title of the plot, default: 'Variance plot' 
  } 
  \item{cex.lab}{ 
    character expansion for column and rownames, default: 1.5 
  } 
  \item{cex.sym}{ 
    character expansion for the plotting symbols, default: 8 
  } 
  \item{cex.var}{ 
    character expansion for the variance values, default: 1 
  } 
  \item{pch}{ 
    plotting character for the variances, default: 16 (filled circle) 
  } 
  \item{...}{ 
    remaining arguments are delegated to the standard plot function 
  } 
} 
\examples{ 
data(USArrests) 
pbi_pca.varplot(prcomp(USArrests,scale=TRUE),cex.sym=8) 
data(swiss) 
pca=prcomp(swiss,scale=TRUE) 
round(pbi_pca.variances(pca),3) 
pbi_pca.varplot(pca,cex.sym=5,cex.var=0.7,cex.lab=1) 
} 
\seealso{  See also [pbi](#pbi), [pbi$pca.variances](#pca.variances) } 
