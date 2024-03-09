\name{pbi_protscale} 
\alias{pbi$protscale} 
\alias{pbi_protscale} 
\title{Calculate protscale moving averages and optionally plot them.} 
\description{ 
  This functions takes a given input sequence and calculates the hydophobicity averages 
  over a slising window of length 9 using the Kyte-Doolitle scale. The averaged scores 
  can be as well plotted. 
} 
\usage{pbi_protscale(sequence,plot=FALSE,col='orange')} 
\arguments{ 
  \item{sequence}{ 
    either a sequence in FASTA format or a simple sequence or text string. 
  } 
  \item{plot}{ 
    should the biophysical properties be plotted, default: FALSE 
  } 
  \item{col}{ 
    color for the plot  
  } 
} 
\value{return vector of window averages (invisible).} 
\examples{ 
fasta=" 
>sp|P04156.1|PRIO_HUMAN 
MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHG 
GGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIH 
FGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVV 
EQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG" 
pbi_protscale(fasta,plot=TRUE,col="light blue") 
} 
