\name{pbi_grect} 
\alias{pbi$grect} 
\alias{pbi_grect} 
\title{Adds colored background and grid lines to an existing plot} 
\description{The function can be used to add a background color to existing plots. 
  In case the plotting was already done you should add a very transparent color,  
  using for instance 33 as the last two digits for a RGB color. 
} 
\usage{ pbi_grect(col='#c0c0c033',grid=FALSE) } 
\arguments{ 
  \item{ col }{ 
    - background color for the plot, default: '#c0c0c033' 
  } 
  \item{ grid }{ 
    - should a grid been drawn, default: FALSE 
  } 
} 
\examples{par(mfrow=c(1,2),mai=rep(0.8,4)) 
  data(iris) 
  pbi_xyplot(iris[,1:2],col=as.numeric(iris$Species)+1) 
  pbi_grect() 
  plot(iris[,1:2],type="n") 
  pbi_grect(col="#ffe0e0",grid=TRUE) 
  points(iris[,1:2],col=as.numeric(iris$Species)+1,pch=19) 
} 
\seealso{  See also: [pbi](#pbi), [pbi_xyplot](#xyplot) } 
