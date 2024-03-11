\name{pbi_report.pval} 
\alias{pbi$report.pval} 
\alias{pbi_report.pval} 
\title{Return a p-value for reporting} 
\description{ 
    Return a p-value for reporting, either giving the three alpha thresholds,  
  <0.05, <0.01, or <0.001 or using the star syntax.  
} 
\usage{pbi_report.pval(p.val,star=FALSE)} 
\arguments{ 
  \item{p.val}{ 
    a numerical p-value. 
  } 
  \item{star}{ 
    boolean, should the one-three star syntax be used, default: FALSE. 
  } 
} 
\value{return pvalue in shown in alpha-threshold or star syntax. If the p-value is not significant, either the value or any empty string is returned if the star syntax is used.} 
\examples{ 
  report.pval = pbi_report.pval 
  report.pval(1/10000) 
  report.pval(1/10000,star=TRUE) 
  report.pval(0.02,star=TRUE) 
  report.pval(0.12,star=TRUE) 
  report.pval(c(0.001,0.01,0.3,0.02)) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_report.chisq.test]{pbi_report.chisq.test}}, \code{\link[pbi:pbi_report.conf.int]{pbi_report.conf.int}} 
} 
