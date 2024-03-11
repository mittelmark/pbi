\name{pbi_dpairs.legend} 
\alias{pbi$dpairs.legend} 
\alias{pbi_dpairs.legend} 
\title{Adding legend top or bottom to a pbi_dpairs or pairs plot} 
\description{ 
    The function \code{pbi_dpairs.legend} allows the user to place a legend outside of a  
  pairs or dpairs plot. 
} 
\usage{pbi_dpairs.legend(labels, col='grey80', pch=15, cex=1)} 
\arguments{ 
  \item{labels}{ 
    txt labels to be plotted 
  } 
  \item{col}{ 
    colors for the plotting characters 
  } 
  \item{pch}{ 
    plotting symbol, default 15 
  } 
  \item{cex}{ 
    the character expansion for text and plotting characters, default: 1 
  } 
} 
\examples{ 
  data(iris) 
  par(omi = c(0.8, 0.4,0.8,0.4)) # reserve some space top and bottom 
  pbi_dpairs(iris,col.box=2:4,col.xy=rep(c(2:4),each=50)) 
  pbi_dpairs.legend(levels(iris$Species),col=2:4) 
  mtext('Iris Data',side=3,outer=TRUE,cex=2,line=1) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_dpairs]{pbi_dpairs}}} 
 
