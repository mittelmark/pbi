\name{pbi_mw} 
\alias{pbi$mw} 
\alias{pbi_mw} 
\title{ Determine the molecular weight for a given sequence (sequences). } 
\description{Calculates and returns the molecular weight for a given amino acid sequence.} 
\usage{ pbi_mw(seq)} 
\arguments{ 
  \item{ seq }{ 
    an amino acid sequence string 
  } 
} 
\value{return the protein molecular weight.} 
\examples{ 
  pbi_mw("AACTLIS") 
  pbi_mw("A") 
} 
