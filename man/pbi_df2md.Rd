\name{pbi_df2md} 
\alias{pbi$df2md} 
\alias{pbi_df2md} 
\title{Convert a data frame or a matrix into a Markdown table} 
\description{ 
  This function can be used within Rmarkdown documents to display easily 
  a simple Markdown table. For more advance use  cases you should other commands 
  such as kable from the knitr package. 
} 
\usage{pbi_df2md(x,caption='',rownames=TRUE)} 
\arguments{ 
  \item{x}{ 
    matrix or data frame  
  } 
  \item{rownames}{ 
    should  the rownames be displayed, default: TRUE 
  } 
  \item{caption}{ 
    the caption for the table, it is just displayed below of the table. 
  } 
} 
\value{return prints to stdout.} 
\examples{ 
  data(swiss) 
  pbi_df2md(head(swiss)) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}} 
 
