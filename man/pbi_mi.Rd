\name{pbi_mi} 
\alias{pbi$mi} 
\alias{pbi_mi} 
\title{ Return the mutual information for two vectors or a binned table. } 
\description{ 
    Returns: mutual information value, or matrix of all pairwise values in case input  is matrix or data frame 
} 
\usage{ pbi_mi(x,y=NULL, breaks=4) } 
\arguments{ 
  \item{ x }{ 
    either a binned table, a numerical vector, a data frame or matrix 
  } 
  \item{ y }{ 
    a numerical vector if x is not a binned table, matrix or data frame 
  } 
  \item{ breaks }{ 
    number of breaks to create a binned table if x and y are numerical vectors, default: 4 
  } 
} 
\value{return mutual information value or matrix of all pairwise values in case input is matrix or data frame.} 
\examples{ 
rn1=rnorm(100,mean=10,sd=1); 
rn2=rn1+0.5*rnorm(100) 
cor(rn1,rn2) # high 
cor(rn1,sample(rn2)) #random  
pbi_mi(rn1,rn2) # high  
pbi_mi(rn1,sample(rn2)) #random 
tab=table(cut(rn1,breaks=4),cut(rn2,breaks=4)) 
pbi_mi(tab) 
pbi_mi(rn1,rn2,breaks=4) 
pbi_mi(rn1,rn2,breaks=7) 
data(iris) 
round(pbi_mi(iris[,1:4]),2) 
mii=pbi_mi(iris[,1:4]) 
round(mii/diag(mii),2) 
round(cor(iris[,1:4]),2) 
} 
