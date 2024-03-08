\name{testprint} 
\alias{testprint} 
\title{ print a test message } 
\description{ 
    A test function which prints a string. 
} 
\usage{ testprint(txt) } 
\arguments{ 
  \item{txt}{ some value to print  } 
} 
\examples{ 
    testprint("Hello World!") 
    # sample file use 
   dec=read.table(file.path(system.file(package="Rpkg"),"files","decathlon.tab")) 
   head(dec) 
} 
\seealso{  See also \code{\link{add}} } 
\keyword{ arith } 
