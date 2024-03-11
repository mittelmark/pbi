\name{pbi_readPepinfo} 
\alias{pbi$readPepinfo} 
\alias{pbi_readPepinfo} 
\title{Read data from the EMBOSS pepinfo tool} 
\description{ 
  This is a function to visualize the biophysical properties of a protein 
  using the output of the EMBOSS pepinfo tool which can be accessed online 
  at \url{https://www.ebi.ac.uk/jdispatcher/seqstats/emboss_pepinfo} 
After submitting your sequence you have to use the file at  
"Result Files->Tool Output" the file ending with ".output". 
} 
\usage{pbi_readPepinfo(file,region=NULL)} 
\arguments{ 
  \item{file}{ 
    the result file from the EMBOSS pepinfo file. 
  } 
  \item{region}{ 
    string matching a certain region such as "Doolittle", if not given all available regions will be shown,  default: NULL 
  } 
} 
\value{return dataframe with the columns: Position, Aminoacid, Result.} 
\examples{ 
pepfile=system.file("files/pepinfo-spike-sars2.txt",package="pbi") 
if (file.exists(pepfile)) { 
  pbi_readPepinfo(pepfile) 
  res=pbi_readPepinfo(pepfile,region="Doolittle") 
  plot(res$Result ~ res$Position,type="h",col="orange") 
} 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_readFasta]{pbi_readFasta}} 
} 
