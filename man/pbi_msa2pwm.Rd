\name{pbi_msa2pwm} 
\alias{pbi$msa2pwm} 
\alias{pbi_msa2pwm} 
\title{Calculates PFM, PPM and PWM for the given alignment file.} 
\description{ 
    Implementation is partially derived from Dave Tang, see: https://davetang.org/muse/2013/10/01/position-weight-matrix/ 
} 
\usage{pbi_msa2pwm(filename)} 
\arguments{ 
  \item{filename}{ 
    filename of a file wither with two columns (second column with aligned sequences, or a single column file with the first column having the sequences. 
  } 
} 
\value{return list with the following components: 
  - PFM - Position frequency matrix 
  - PPM - Position probability matrix 
  - PWM - Position weight matrix 
  - PWMPC - Position weight matrix with added pseudocounts (untested)} 
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
