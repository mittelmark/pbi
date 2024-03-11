\name{pbi_dpairs} 
\alias{pbi$dpairs} 
\alias{pbi_dpairs} 
\title{Improved pairs plot considering the data types} 
\description{ 
  The function \code{pbi_dpairs} provides an improved pairs plot which accounts 
  for the data type of the actual variables. It will plot in the  
  lower diagonal xy-plots, box-plots or assoc-plots depending on the  
  two data types. In the upper diagonal effect sizes and stars for the p-values 
  for the tests (anova, t.test, chisq.test or cor.test) will be shown. In the diagonal  
  the data distribution will be outlined. This plot is usually an useful visualization for 3-8 variables. 
} 
\usage{pbi_dpairs(data, col.box='grey80', col.xy="grey60", cex.diag=2.5, order=TRUE,  
    pch=19)} 
\arguments{ 
  \item{data}{ 
    data frame with columns of class factor, numeric or integer. 
  } 
  \item{col.box}{ 
    colors for the boxplots, either a single value or a vector of colors for each level of a factor variable, default; 'grey80' 
  } 
  \item{col.xy}{ 
    colors for the xy-plots, either a single value of a vector which is as long as the number of data points, default: 'grey60' 
  } 
  \item{cex.diag}{ 
    character expansion for the diagonal texts 
  } 
  \item{order}{ 
    should the variables be ordered by data type and name, this is recommended as it orders the plots, starting with assocplots, then boxplots and finally xyplots, default: TRUE 
  } 
  \item{pch}{ 
    plotting character for xy-plots, default 19 (round circle). 
  } 
} 
\examples{ 
  data(iris) 
  par(omi = c(0.8, 0.4,0.4,0.4)) 
  pbi_dpairs(iris,col.box=2:4,col.xy=rep(c(2:4),each=50), 
     cex.diag=1.6) 
  pbi_dpairs.legend(levels(iris$Species),col=2:4) 
 
  par(omi=c(0.5,0.5,0.8,0.2)) 
  library(MASS) 
  btwt=birthwt;  
  for (col in c('low','race','smoke','ptl','ht','ui','ftv')) {  
     btwt[,col]=as.factor(btwt[,col])  
  } 
  pbi_dpairs(btwt[,2:8],cex.diag=1.6) 
      mtext('Birth-Weight data',side=3,outer=TRUE, 
      cex=1.5,line=1) 
  if (require("palmerpenguins")) { 
      data(penguins) 
      peng=penguins 
      colnames(peng)[3]='bill.len' 
      colnames(peng)[4]='bill.dep' 
      colnames(peng)[5]='flip.len' 
      colnames(peng)[6]='mass' 
      pbi_dpairs(peng,col.xy=pbi_pastel(3)[as.numeric(peng$species)]) 
      pbi_dpairs.legend(levels(peng$species), 
          col=pbi_pastel(3)[as.numeric(as.factor(levels(peng$species)))]) 
  } 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_dpairs.legend]{pbi_dpairs.legend}}} 
 
