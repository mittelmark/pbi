\name{pbi_wordFreq} 
\alias{pbi$wordFreq} 
\alias{pbi_wordFreq} 
\title{Count the number of words with a given length in a sequence.} 
\description{ 
  The function creates a sliding window of length `wlength` over the given  
  text or sequence string creating words of length `wlength`.  
  The number how often a certain word appears is counted.  
  So - `AABAAC` - contains the words: `AA`, `AB`, `BA`, `AA` and `AC`. 
} 
\usage{pbi_wordFreq(seq,wlength=2)} 
\arguments{ 
  \item{seq}{ 
    sequence as text string or vector of text strings 
  } 
  \item{wlength}{ 
    word size, default: 2 
  } 
} 
\value{return list object with words as keys and counts as values or data frame in case input is a vector.} 
\examples{ 
  seq="AAABBBCCCDEFAABBCCDD" 
  unlist(pbi_wordFreq(seq)) 
  unlist(pbi_wordFreq(seq,wlength=4)) 
  minions=read.table(text=' 
  Minion1      STSTTS 
  Minion2      TTTTTT 
  Minion3      STSTTT 
  Minion4      STTTTT 
  Minion5      SSTTTT 
  Minion6      SSSTST 
  Minion7      SSSSTT 
',row.names=1) 
  min=pbi_wordFreq(minions[,1],wlength=2) 
  rownames(min)=rownames(minions) 
  min 
} 
