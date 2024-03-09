\name{pbi_pca.corplot} 
\alias{pbi_pca.corplot} 
\alias{pbi$pca.corplot} 
\title{PCA correlation plot to show association between PCs and variables.} 
\description{ 
  The function provides a PCA correlation plot to show associations  
  between PCs and variables. The closer a variable to the PC coordinate  
  the higher the correlation, the more away from the center of the  
  coordinate system, the higher the impact of the variable on this PC. 
  You can think about the corplot as a biplot ommiting the samples and  
  the arrows. 
} 
\usage{pbi_pca.corplot(pca,pcs=c("PC1","PC2"), main="Correlation plot",cex=NULL,nvar=64,...)} 
\arguments{ 
  \item{pca}{ 
    pca object which was created using the function _prcomp_. 
  } 
  \item{pcs}{ 
    vector of two PCs to be plotted against each other, default: c('PC1','PC2') 
  } 
  \item{main}{ 
    title of the plot, default: 'Correlation plot' 
  } 
  \item{cex}{ 
    character expansion for the samples, default: NULL (automatic calculation) 
  } 
  \item{nvar}{ 
    number of variables which will be displayed, for both components the variables which contributes mostly to the variances will be used, default: 64. 
  } 
  \item{...}{ 
    remaining arguments are delegated to the standard plot function 
  } 
} 
\examples{ 
par(mfrow=c(1,2),mai=c(0.7,0.7,0.5,0.1)) 
library(cluster) 
data(votes.repub) 
pca=prcomp(t(na.omit(votes.repub))) 
pbi_pca.corplot(pca) 
data(swiss) 
pca=prcomp(swiss,scale.=TRUE) 
pbi_pca.corplot(pca) 
} 
