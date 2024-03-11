\name{pbi_msa2pwm} 
\alias{pbi$msa2pwm} 
\alias{pbi_msa2pwm} 
\title{PFM, PPM and PWM for a alignment file} 
\description{ 
    Implementation is partially derived from Dave Tang, see: \url{https://davetang.org/muse/2013/10/01/position-weight-matrix/} 
} 
\usage{pbi_msa2pwm(filename)} 
\arguments{ 
  \item{filename}{ 
    filename of a file wither with two columns (second column with aligned sequences, or a single column file with the first column having the sequences. 
  } 
} 
\value{return list with the following components: 
\itemize{ 
  \item{PFM - Position frequency matrix} 
  \item{PPM - Position probability matrix} 
  \item{PWM - Position weight matrix} 
  \item{PWMPC -  Position weight matrix with added pseudocounts (untested)} 
} 
} 
\examples{ 
# data example: https://en.wikipedia.org/wiki/Position_weight_matrix 
  cat("GAGGTAAAC 
TCCGTAAGT 
CAGGTTGGA 
ACAGTCAGT 
TAGGTCATT 
TAGGTACTG 
ATGGTAACT 
CAGGTATAC 
TGTGTGAGT 
AAGGTAAGT 
",file="test.clu") 
lapply(pbi_msa2pwm("test.clu"),round,2) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_mw]{pbi_mw}}, \code{\link[pbi:pbi_readFasta]{pbi_readFasta}} 
} 
