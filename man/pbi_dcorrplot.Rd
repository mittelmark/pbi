\name{pbi_dcorrplot} 
\alias{pbi$dcorrplot} 
\alias{pbi_dcorrplot} 
\title{Visualize a correlation matrix.} 
\description{ 
    Visualize a correlation matrix. 
} 
\usage{pbi_dcorrplot(mt,text.lower=TRUE, text.upper=FALSE,pch=19,p.mat=NULL,alpha=0.05,cex.sym=5,cex.r=1,cex.lab=1.4,...)} 
\arguments{ 
  \item{mt}{ 
    matrix with pairwise correlations 
  } 
  \item{text.lower}{ 
    should in the lower diagonal the correlation coefficient be shown, default: TRUE 
  } 
  \item{text.upper}{ 
    should in the upper diagonal the correlation coefficient be shown, default: FALSE 
  } 
  \item{pch}{ 
    the plotting symbol for the correlations, default: 19 
  } 
  \item{p.mat}{ 
    matrix with p-values to strike out insignificant p-values, default: NULL (not used) 
  } 
  \item{alpha}{ 
    significance threshold for _p.mat_, default: 0.05 
  } 
  \item{cex.sym}{ 
    character expansion for the correlation symbols, default: 5 
  } 
  \item{cex.r}{ 
    character expansion for the r-values if _text.lower_ or _text.upper_ are set to TRUE, default: 1 
  } 
  \item{cex.lab}{ 
    character expansion for the variable text labels, default: 1.4 
  } 
  \item{...}{ 
    arguments delegated to the plot function 
  } 
} 
\examples{ 
  data(swiss) 
  sw=swiss 
  colnames(sw)=abbreviate(colnames(swiss),6) 
  cr=pbi_dcorr(sw,method='spearman') 
  pbi_dcorrplot(cr$estimate,cex.sym=8,text.lower=TRUE, 
     cex.r=1.5,p.mat=cr$p.value) 
} 
\seealso{  See also [pbi](#home), [pbi_dcorr](#dcorr) } 
