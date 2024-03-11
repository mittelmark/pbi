\name{pbi_cohensH} 
\alias{pbi$cohensH} 
\alias{pbi_cohensH} 
\title{Effect size for a 2x2 contingency table} 
\description{ 
  The function \code{pbi_cohensH} calculates the effect size for 2x2 contingency tables.  
  Due to Cohen's rule of thumb values of 0.2 to 0.5 are considered to stand  
  for small effects, values from 0.5 to 0.8 represent medium effects and values above 0.8 represent large effects.  
} 
\usage{pbi_cohensH(tab)} 
\arguments{ 
  \item{tab}{ 
    a 2x2 contingency table 
  } 
} 
\value{return Cohen's h value} 
\examples{ 
# data from New Eng. J. Med. 329:297-303, 1993 
azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2)) 
rownames(azt)=c("AZT","Placebo") 
colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress") 
pbi_cohensH(azt) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_cohensW]{pbi_cohensW}} } 
