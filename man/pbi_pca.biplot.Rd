\name{pbi_pca.biplot} 
\alias{pbi$pca.biplot} 
\alias{pbi_pca.biplot} 
\title{Improved biplot for prcomp objects} 
\description{ 
  The function \code{pca.biplot} provides an improved biblot for 
  visualizing the pairwise scores of individual principal components of  
  an object created using the function _prcomp_. In contrast to the default  
  biplot function  this plot visualizes the data as points and not row numbers, 
  it allows to display groups using color codes and distribution ellipses. 
} 
\usage{pbi_pca.biplot(pca,pcs=c("PC1","PC2"), pch=19, col='black', text=NULL, 
    arrows=TRUE, arrow.fac=1, arrow.n=-1, ellipse=FALSE, ell.fill=FALSE,  
    xlab=NULL, ylab=NULL, grid=TRUE, scale=NULL, ...)} 
\arguments{ 
  \item{pca}{ 
    pca object of class _prcomp_, created using the function _prcomp_. 
  } 
  \item{pcs}{ 
    the components to plot, default: c('PC1','PC2') 
  } 
  \item{pch}{ 
    plotting character, default: 19 
  } 
  \item{col}{ 
    plotting color, default: black 
  } 
  \item{text}{ 
    optional text labels, if not given  pch is used, default: NULL 
  } 
  \item{arrows}{ 
    should loading arrows be displayed, default: TRUE 
  } 
  \item{arrow.fac}{ 
    scaling factor for arrow length, default: 1 
  } 
  \item{arrow.n}{ 
    how many arows per PC to display, default: -1 (all) 
  } 
  \item{ellipse}{ 
    should 85 and 95 confidence intervals for the chisq distribution be shown. If this is shown colors for each group using the col argument must be given, default: FALSE 
  } 
  \item{ell.fill}{ 
    should a filled 85 percent confidence interval be shown, colors will be used from the plotting color with opacity, default: FALSE 
  } 
  \item{xlab}{ 
    custom xlab, if not given the PC name with variance in percent is shown, default: NULL 
  } 
  \item{ylab}{ 
    custom ylab, if not given the PC name with variance in percent is shown, default: NULL 
  } 
  \item{grid}{ 
    should a plot grid be drawn, default: TRUE 
  } 
  \item{scale}{ 
    function to scale to coordinate values sich as asinh, default: NULL 
  } 
  \item{...}{ 
    additional arguments delegated to the standard plot function 
  } 
} 
\examples{ 
  par(mfrow=c(1,2),mai=c(0.8,0.8,0.6,0.2)) 
  data(iris) 
  pci=prcomp(iris[,1:4],scale=TRUE) 
  pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=TRUE,ell.fill=TRUE, 
      arrow.fac=2.3,arrows=TRUE) 
  legend('topright',pch=19,col=2:4,levels(iris$Species)) 
  # just a score plot 
  pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=TRUE,ell.fill=TRUE, 
      arrow.fac=2.3,arrows=FALSE) 
  pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=FALSE,arrow.fac=2.3,arrows=FALSE) 
  pbi_pca.biplot(pci,pcs=c('PC1','PC3'),col=rep(2:4,each=50), 
         ellipse=FALSE,arrow.fac=2.3,arrows=FALSE) 
  data(swiss) 
  col=c(2,4)[as.numeric(cut(swiss$Catholic,breaks=c(0,20,100)))] 
  pcs=prcomp(swiss,scale=TRUE) 
  pbi_pca.biplot(pcs,arrow.fac=2,grid=FALSE, 
        col=col,ellipse=TRUE,text=substr(rownames(swiss),1,3)) 
  pbi_pca.biplot(pcs,arrow.fac=2,grid=TRUE,  
       col=col,text=substr(rownames(swiss),1,3),scale=asinh,arrows=FALSE) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_pca.corplot]{pbi_pca.corplot}},  
\code{\link[pbi:pbi_pca.pairs]{pbi_pca.pairs}},  
\code{\link[pbi:pbi_pca.plot]{pbi_pca.plot}}, 
\code{\link[pbi:pbi_pca.toData]{pbi_pca.toData}}, 
\code{\link[pbi:pbi_pca.variances]{pbi_pca.variances}}, 
\code{\link[pbi:pbi_pca.varplot]{pbi_pca.varplot}} 
} 
