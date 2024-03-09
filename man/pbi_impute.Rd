\name{pbi_impute} 
\alias{pbi$impute} 
\alias{pbi_impute} 
\title{ Missing value imputation using mean, rpart or knn methods. 
} 
\description{ 
    Replaces missing values with a reasonable guess by different imputation methods such as  
  the simple and not recommended methods mean and median, where NA's are replaced with the  
  mean or median for this variable or the more recommended methods using rpart decision trees 
  or knn using a correlation distance based k-nearest neighbor approach.  
  The rpart method can be as well used to replace missing values for categorical variables. 
  In case of median and mean imputations for categorical variables the modus is used,  
  so missing values are replaced with the most often category. This is rarely reasonable. 
} 
\usage{pbi_impute(x,method="mean", k=5, cor.method="spearman")} 
\arguments{ 
  \item{x}{either a matrix or data frame} 
  \item{method}{ the method used for replacing missing values, either mean,  
     median, rpart or knn, default mean (not recommended). 
  } 
  \item{k}{in case of knn imputation the number of neighbors, default: 5} 
  \item{cor.method}{in case of distance determination the correlation method, default: 'spearman'} 
} 
\value{return depending on the input either a data frame or matrix with NA's replaced by imputed values.} 
\examples{ 
  data(iris) 
  ir=as.matrix(iris[,1:4]) 
  idx=sample(1:length(ir),50) 
  ir[idx]=NA 
  summary(ir) 
  irc=as.matrix(ir) 
  ir=pbi_impute(ir) 
  summary(ir) 
  ir=iris 
  ir[1,3]=NA; ir[2,4]=NA; ir[c(1,3),5]=NA 
  head(ir,3) 
  head(pbi_impute(ir,method="rpart"),3) 
  irmea=pbi_impute(irc,method="mean") 
  irknn=pbi_impute(irc,method="knn") 
  irrpt=pbi_impute(irc,method="rpart") 
  cor(as.matrix(iris[,1:4])[is.na(irc)],irmea[is.na(irc)])   
  cor(as.matrix(iris[,1:4])[is.na(irc)],irknn[is.na(irc)])   
  cor(as.matrix(iris[,1:4])[is.na(irc)],irrpt[is.na(irc)])   
} 
\seealso{  See also: [pbi](#pbi) } 
