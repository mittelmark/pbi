\name{pbi_report.conf.int} 
\alias{pbi$report.conf.int} 
\alias{pbi_report.conf.int} 
\title{Formatted text strings for reporting a confidence interval} 
\description{ 
    Return a formatted text string for reporting a confidence interval.  
} 
\usage{pbi_report.conf.int(ci,round=2)} 
\arguments{ 
  \item{ci}{ 
    a confidence interval consisting of two numerical values 
  } 
  \item{round}{ 
    rounding digits, default: 2 
  } 
} 
\value{return a formatted text string for reporting a confidence interval to be included in a LaTeX/Sweave document.} 
\examples{ 
  azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2)) 
  rownames(azt)=c("AZT","Placebo") 
  colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress") 
  pbi_report.conf.int(prop.test(azt)$conf.int) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_report.chisq.test]{pbi_report.chisq.test}}, \code{\link[pbi:pbi_report.pval]{pbi_report.pval}} 
} 
