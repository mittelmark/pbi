\name{pbi_file.head} 
\alias{pbi$file.head} 
\alias{pbi_file.head} 
\title{Displays the first n lines of a file to the terminal} 
\description{ 
    Displays the first n lines of a file to the terminal 
} 
\usage{pbi_file.head(filename,n=6)} 
\arguments{ 
  \item{filename}{ 
    filename of a text file 
  } 
  \item{n}{ 
    number of first lines to display, default: 6 
  } 
} 
\examples{ 
  pbi_file.head(system.file( 
  "files/pepinfo-spike-sars2.txt",package="pbi")) 
} 
\seealso{  See also [sbi$file.cat](#file.cat) } 
