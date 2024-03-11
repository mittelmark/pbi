\name{pbi_dassoc} 
\alias{pbi$dassoc} 
\alias{pbi_dassoc} 
\title{Assocplots with residual coloring} 
\description{ 
   This function updates the standard assocplot function from the graphics package  
   with the ability to display residual colors. In blue and red are shown groups with  
   residuals above +4 or below -4 in light colors are shown residuals between 2 and 4 for positive and -4 and -2 for negative residuals. 
} 
\usage{pbi_dassoc(...,shade=TRUE)} 
\arguments{ 
  \item{...}{ 
    arguments delegated to the standard assocplot function. 
  } 
  \item{shade}{ 
    should the residuals been shown, default: TRUE 
  } 
} 
\examples{ 
  x = margin.table(HairEyeColor, c(1, 2)) 
  pbi_dassoc(x) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}} 
 
