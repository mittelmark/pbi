\name{pbi_dist} 
\alias{pbi$dist} 
\alias{pbi_dist} 
\title{Distance function which as well supports binary and correlation distances in addition to standard distances.} 
\description{ 
  This function is an extension to stats::dist function as  it supports as well correlation distance and  
  binary distance measures such as Jaccard coefficient and Matching coefficient.  
  The correlation distance is implemented as `D(i,j) = 1-((r(i,j)+1)/2)`  
  so negative correlations have low similarities. 
} 
\usage{pbi_dist(x,method="euclidean",...)} 
\arguments{ 
  \item{x}{ 
    data frame or matrix with numerical data,  
    in case of binary data as well boolean and two level nominal data could be supplied. 
  } 
  \item{method}{ 
    the distance measure to be used, one of the measures  
    for `stats::dist` such as "euclidean" or "correlation",  
    alias for  "pearson" or "spearman","kendall",  
    for binary data "jc" (Jaccard) and "mc" (Matching coeffient) are supported. 
  } 
  \item{...}{ 
    remaining arguments are forwarded to stats::dist function in case the method is  
    handled by this default method.  
  } 
} 
\value{return distance matrix} 
\examples{ 
  round(cor(iris[,1:4]),2) 
  round(as.matrix(pbi_dist(t(iris[,1:4]),method="pearson")),2) 
  biris=iris[,1:4] 
  biris=apply(biris,2,function(x) { return(x>median(x)) }) 
  head(biris,3) 
  head(apply(biris,2,as.numeric),3) 
  biris=apply(biris,2,as.numeric) 
  summary(biris) 
  round(as.matrix(pbi_dist(t(biris),method="mc")),2) 
  d.can=as.matrix(pbi_dist(t(scale(iris[,1:4])),method="can")) 
  d.can 
  d.can=d.can/max(d.can) 
  round(d.can,2) 
} 
