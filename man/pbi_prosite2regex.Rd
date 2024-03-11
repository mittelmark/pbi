\name{pbi_prosite2regex} 
\alias{pbi$prosite2regex} 
\alias{pbi_prosite2regex} 
\title{Convert PROSITE patterns to regular expressions} 
\description{ 
    This little function converts a PROSITE pattern to a normal regular  
  expression. Please note, that for searching in FASTA files those  
  regular expressions did work correctly if the sequence pattern 
  in the FASTA file is splitted over two files, the sequences should be 
  before fused to single line sequences. The function  
  [pbi_searchFasta](#searchFasta) does this. 
} 
\usage{pbi_prosite2regex(pattern)} 
\arguments{ 
  \item{pattern}{ 
    a PROSITE pattern. 
  } 
} 
\value{return a "normal" regular expression.} 
\examples{ 
  pbi_prosite2regex("A-T-x(0,3)-{ALV}-A") 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}} 
} 
