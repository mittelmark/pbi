\name{pbi_cohensW} 
\alias{pbi$cohensW} 
\alias{pbi_cohensW} 
\title{Effect size for contingency tables} 
\description{ 
  The function \code{pbi_cohensW} calculates the effect size for contingency tables.  
  Due to Cohen's rule of thumb values of 0.1 to 0.3 are considered to stand  
  for small effects, values from 0.3 to 0.5 represent medium effects and values  
  above 0.5 represent large effects. 
} 
\usage{pbi_cohensW(x,p=NULL)} 
\arguments{ 
  \item{x}{ 
    contingency table with counts, usually created using the table  
    command for two variables  or vector with counts with observations 
    for a single variable, if x is a vector p must be given 
  } 
  \item{p}{ 
    expected proportions if x is a vector, in case *x* has length of two, 
    a single value can be given, otherwise p must have same length as x, default: NULL 
  } 
} 
\examples{ 
data(Titanic) 
Titanic[1,1,,] 
pbi_cohensW(Titanic[1,1,,]) 
# Data from New Eng. J. Med. 329:297-303, 1993 
azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2)) 
rownames(azt)=c("AZT","Placebo") 
colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress") 
pbi_cohensW(azt) 
# number of boys (25) and girls (15) in a hospital which deviates  
# from 50/50 or 51/49 ratios 
pbi_cohensW(c(25,15),p=0.5) 
pbi_cohensW(c(25,15),p=c(0.51,0.49)) 
# most extrem case 40 boys and 0 girls 
pbi_cohensW(c(40,0),p=0.5) 
pbi_cohensW(c(40,0),p=c(0.51,0.49)) # max value here around 2*0.49 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_cohensH]{pbi_cohensH}} } 
