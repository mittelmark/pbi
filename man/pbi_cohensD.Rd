\name{pbi_cohensD} 
\alias{pbi$cohensD} 
\alias{pbi_cohensD} 
\title{Effect size for the difference between two means} 
\description{ 
  The function \code{pbi_cohensD} calculates the effect size for the difference between two means. 
  Due to Cohen's rule of thumb values of 0.2 to 0.5 are considered to stand  
  for small effects, values from 0.5 to 0.8 represent medium effects and values above 0.8 represent large effects. 
} 
\usage{pbi_cohensD(x,y,paired=FALSE)} 
\arguments{ 
  \item{x}{ 
    vector with numerical values 
  } 
  \item{y}{ 
    vector with two grouping variables, having the same length as num 
  } 
  \item{paired}{ 
    are the data paired, default: FALSE 
  } 
} 
\value{return Cohen's d value} 
\examples{ 
  set.seed(125) 
  data(sleep) 
  with(sleep,pbi_cohensD(extra,group)) 
  x1=rnorm(100,mean=20,sd=1) 
  x2=rnorm(100,mean=22,sd=1) 
  g1=rep('A',100) 
  g2=rep('B',100) 
  # difference should be around 2SD 
  pbi_cohensD(c(x1,x2),as.factor(c(g1,g2))) 
  # biserial correlation coefficient as alternative 
  # value is as well large 
  cor(c(x1,x2),as.numeric(as.factor(c(g1,g2)))) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, \code{\link[pbi:pbi_cohensW]{pbi_cohensW}}  
\code{\link[pbi:pbi_etaSquared]{pbi_etaSquared}}, \code{\link[pbi:pbi_epsilonSquared]{pbi_epsilonSquared}} } 
