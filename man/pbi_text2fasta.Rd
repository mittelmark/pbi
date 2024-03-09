\name{pbi_text2fasta} 
\alias{pbi$text2fasta} 
\alias{pbi_text2fasta} 
\title{Convert any Text file to a FASTA file using only Aminoacid letters.} 
\description{ 
    This function loops over the given directory and converts text files to 
  a FASTA file format. Only uppercase letters are used for the output file. 
} 
\usage{pbi_text2fasta(dir,pattern="*",outfile="stdout")} 
\arguments{ 
  \item{dir}{ 
    a directory 
  } 
  \item{pattern}{ 
    regular expression for the files which should be analysed 
  } 
  \item{outfile}{ 
    the output filename, default: 'stdout' 
  } 
} 
\value{return the number of sequences.} 
\examples{ 
  pbi_text2fasta(dir=".",pattern='*.R$',outfile="test.fasta") 
  pbi_file.head("test.fasta",n=3) 
} 
\seealso{  See also [pbi](#home), [pbi_searchFasta](#searchFasta) } 
