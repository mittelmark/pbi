\name{pbi_dcorr} 
\alias{pbi$dcorr} 
\alias{pbi_dcorr} 
\title{Pairwise correlations and their p-values} 
\description{ 
  The function is an extension to the standard \code{\link[stats:cor]{stats::cor}} function, it calculates as well 
  the p-values for the pairwise associations and returns them in a matrix as well. 
} 
\usage{pbi_dcorr(data,method='pearson',use='pairwise.complete.ob')} 
\arguments{ 
  \item{data}{ 
    matrix or data frame where the variables are in the columns, NA's are allowed. 
  } 
  \item{method}{ 
    type of correlation to be determined, either 'pearson', 'spearman' or 'kendall', default: 'pearson' 
  } 
  \item{use}{ 
    how to deal with NA's, default: 'pairwise.complete.obs' 
  } 
} 
\value{returns list with the following components: 
  \itemize{ 
     \item \emph{estimate} - matrix with correlation values 
     \item \emph{p.value}  - matrix with p-values 
     \item \emph{method} - character string with the used correlation method 
  } 
} 
\examples{ 
  data(swiss) 
  lapply(pbi_dcorr(swiss)[1:2],round,2) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_dcorrplot]{pbi_dcorrplot}}, \code{\link[pbi:pbi_mi]{pbi_mi}}  
} 
 
