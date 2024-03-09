\name{pbi_epsilonSquared} 
\alias{pbi$epsilonSquared} 
\alias{pbi_epsilonSquared} 
\title{Calculate the effect size epsilon-squared for variables of a kruskal.test} 
\description{ 
    Cohen's rule of thumb for interpretation is: 0.01-0.09 small, 0.09-0.25 medium and above 0.25 large effect. You can convert epsilon-squared to a Spearman _r_ by using the square root of epsilon-square. 
} 
\usage{pbi_epsilonSquared(x,y)} 
\arguments{ 
  \item{x}{ 
    vector with numerical values or a linear model or an aov object. 
  } 
  \item{y}{ 
    vector with factor values or numerical vector of a second class 
  } 
} 
\value{return numerical value for epsilon-square for given variables.} 
\examples{ 
  data(iris) 
  epsilonSquared=pbi_epsilonSquared 
  epsilonSquared(iris$Sepal.Length,iris$Species) 
  # two factor example as well for wilcox.test possible 
  data(ToothGrowth) 
  epsilonSquared(ToothGrowth$len,as.factor(ToothGrowth$dose)) 
  # close to r-square of spearman! 
  cor(ToothGrowth$len,ToothGrowth$dose,method="spearman")^2 
} 
\seealso{  See also [pbi](#home), [pbi_etaSquared](#etaSquared) } 
