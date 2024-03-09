\name{pbi_report.chisq.test} 
\alias{pbi$report.chisq.test} 
\alias{pbi_report.chisq.test} 
\title{Return a formatted text string for reporting a chisq.test.} 
\description{ 
    Return a formatted text string for reporting a chisq.test. 
} 
\usage{pbi_report.chisq.test(tab)} 
\arguments{ 
  \item{tab}{ 
    a contigency table. 
  } 
} 
\value{return formatted text string for reporting a chisq.test in a LaTeX/Sweave document.} 
\examples{ 
  azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2)) 
  rownames(azt)=c("AZT","Placebo") 
  colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress") 
  pbi_report.chisq.test(azt) 
} 
\seealso{  See also [pbi](#home), [pbi_report.pval](#report.pval) } 
