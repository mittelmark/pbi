\name{pbi_wininstall} 
\alias{pbi$wininstall} 
\alias{pbi_wininstall} 
\title{Make an executable Batch script on Windows for R applications within the users PATH.} 
\description{ 
  The function two files in the users PATH for executables on Windows, a BATCH file and a script file 
  both in the same folder and with the same file prefix. This allows the user to directly execute a Rscript  
  by pressing Win-R as a shortcut and the then entering the application name assume you have a file: _hw.R_ 
} 
\usage{pbi_wininstall(filename="")} 
\arguments{ 
  \item{filename}{ 
    optional filename which should be installed directly without asking for the filename, default: "" 
  } 
} 
\examples{ 
  \dontrun{ 
  pbi_wininstall("hw.R") 
} 
} 
