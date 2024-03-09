\name{pbi_searchFasta} 
\alias{pbi$searchFasta} 
\alias{pbi_searchFasta} 
\title{Search a FASTA file with a regular expression.} 
\description{ 
    This function searches FASTA files by removing line breaks within the 
  sequence belonging to the same ID. 
} 
\usage{pbi_searchFasta(filename,pattern)} 
\arguments{ 
  \item{filename}{ 
    a sequence file in FASTA format 
  } 
  \item{pattern}{ 
    a standard regular expression which will be applied on the sequence for the ID. 
  } 
} 
\value{return the matching ID's.} 
\examples{ 
  # the pattern splits over two lines for the first sequence 
  pbi_searchFasta(filename=system.file("files/human-tRNAs.fasta",package="pbi"), 
     pattern="GACGC{5}AT{2}CTCT") 
} 
\seealso{  See also [pbi](#home), [pbi_tkregex](#tkregex), [pbi_prosite2regex](#prosite2regex) } 
