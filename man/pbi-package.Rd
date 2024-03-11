\name{pbi-package} 
\alias{pbi-package} 
\title{The pbi R Package for the Course Practical Bioinformatics} 
\description{The pbi package contains functions used in the course 
Practical Bioinformatics at the University of Potsdam.} 
\details{Some more details: 
The following object is available: 
\describe{ 
\item{\link[pbi:pbi-class]{pbi-class}}{The pbi environment} 
} 
All functions can be called in a method like style using the pbi environment  
or in a function style using the \code{pbi_} - prefix. For example the function \code{pbi_cohensD} can be called as well as \code{pbi$cohensD}. 
} 
\examples{ 
library(pbi) 
ls(pbi) 
ls('package:pbi') 
} 
\author{Detlef Groth <dgroth@uni-potsdam.de>} 
