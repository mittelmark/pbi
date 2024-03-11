\name{pbi_readFasta} 
\alias{pbi$readFasta} 
\alias{pbi_readFasta} 
\title{Read in a (small) fasta file into a list} 
\description{ 
    This function can be used to read in a small FASTA file into a listg where the keys are the ids and the values are 
  are the sequences without line breaks. This code is very slow for large files. 
} 
\usage{pbi_readFasta(filename)} 
\arguments{ 
  \item{filename}{ 
    a sequence file in FASTA format 
  } 
} 
\value{return a list with the sequence ids as keys and the sequences as values.} 
\examples{ 
fout = file('minions.fasta','w') 
cat(">Minion1\nSTSTTS\n>Minion2\nTTTTTT\n>Minion3\nSTSTTT\n",file=fout) 
cat(">Minion4\nSTTTTT\n>Minion5\nSSTTTT\n>Minion6\nSSSTST\n",file=fout) 
cat("\n>Minion7\nSSSSTT\n", file=fout) 
close(fout) 
seq=pbi_readFasta("minions.fasta") 
seq[7] 
M=adist(unlist(seq)) 
M 
plot(hclust(as.dist(M))) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_searchFasta]{pbi_searchFasta}} 
} 
