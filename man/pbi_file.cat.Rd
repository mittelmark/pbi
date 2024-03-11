\name{pbi_file.cat} 
\alias{pbi$file.cat} 
\alias{pbi_file.cat} 
\title{Displays a file to the terminal} 
\description{ 
    Displays the a file to the terminal. Works in a platform independent way. 
} 
\usage{pbi_file.cat(filename)} 
\arguments{ 
  \item{filename}{ 
    filename of a text file 
  } 
} 
\examples{ 
  head(pbi_file.cat(system.file( 
  "files/pepinfo-spike-sars2.txt",package="pbi")),n=10) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_file.head]{pbi_file.head}} 
} 
 
