\name{pbi_pastel} 
\alias{pbi$pastel} 
\alias{pbi_pastel} 
\title{Up to 20 pastel colors} 
\description{ 
  This is an alternative color creation function for R versions before 3.6 where  
  the function `hcl.colors` is not available. 
} 
\usage{ pbi_pastel(n) } 
\arguments{ 
  \item{ n }{ 
    number of colors requested, must be within 2 and 20 
  } 
} 
\value{return vector of colors in RGB codes of requested length 'n'.} 
\examples{ 
pbi_pastel(4) 
par(mai=c(0.2,0.2,0.2,0.1)) 
plot(1:20,col=pbi_pastel(20),cex=3,pch=15) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}} 
} 
