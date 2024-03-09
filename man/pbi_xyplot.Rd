\name{pbi_xyplot} 
\alias{pbi$xyplot} 
\alias{pbi_xyplot} 
\title{Improved xyplot which as well displays a grid and the correlation coefficient.} 
\description{ 
  This is just a simple illustrative function to demonstrate how a standard 
  R function can be modified and extended using different type of arguments 
  overwriting default settings and using the ellipsis and delegating  
  remaining arguments to the default function 
} 
\usage{pbi_xyplot(x,y,col="blue",pch=19,grid=TRUE,xlab=NULL,ylab=NULL,ellipse=FALSE,ell.fill=FALSE,show.r=TRUE,...)} 
\arguments{ 
  \item{x}{ 
    vector with  numerical values, or a matrix or data frame 
  } 
  \item{y}{ 
    vector with  numerical values, can be NULL if x is matrix or data frame, default: NULL 
  } 
  \item{col}{ 
    color for the plotting symbols 
  } 
  \item{pch}{ 
    plotting character, default: 19 (filled circle) 
  } 
  \item{grid}{ 
    should be a grid drawn, default: TRUE 
  } 
  \item{xlab}{ 
    xlabel, default: NULL (x) 
  } 
  \item{ylab}{ 
    ylabel, default: NULL (y) 
  } 
  \item{ellipse}{ 
    draw an ellipse for 85% and 95% confidence intervals using the given colors, default: FALSE 
  } 
  \item{ell.fill}{ 
    should the ellipse be filled for the 85% confidence intervall,  default: FALSE 
  } 
  \item{show.r}{ 
    should the Pearson correlation been shown on top, default: TRUE 
  } 
  \item{...}{ 
    remaining arguments delegated to the standard plot function 
  } 
} 
\examples{ 
data(iris) 
par(mfrow=c(2,2),mai=c(0.8,0.8,0.8,0.1)) 
pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length, 
    xlab="Sepal.Width",ylab="Sepal.Height", 
    col=as.numeric(iris$Species)+1) 
pbi_xyplot(iris$Petal.Width,iris$Petal.Length, 
    xlab="Petal.Width",ylab="Petal.Height", 
    col=as.numeric(iris$Species)+1) 
legend("bottomright",fill=2:4,legend=levels(iris$Species)) 
pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length, 
    xlab="Sepal.Width",ylab="Sepal.Height", 
    col=as.numeric(iris$Species)+1,ellipse=TRUE) 
pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length, 
    xlab="Sepal.Width",ylab="Sepal.Height", 
    col=as.numeric(iris$Species)+1,ellipse=TRUE,ell.fill=TRUE) 
} 
