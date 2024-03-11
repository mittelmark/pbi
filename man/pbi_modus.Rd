\name{pbi_modus} 
\alias{pbi$modus} 
\alias{pbi_modus} 
\title{Most often level in a categorical variable} 
\description{ Simple implementation of returning the most most often appearinglevel(s) of a categorical variable. } 
\usage{ pbi_modus(x) } 
\arguments{ 
  \item{ x }{ 
    a vector with elements of class factor 
  } 
} 
\value{return most often apparent level in the categorical variable} 
\examples{ 
  pbi_modus(c('A','A','B','C')) 
  pbi_modus(c('A','A','B','B','C')) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}} 
} 
