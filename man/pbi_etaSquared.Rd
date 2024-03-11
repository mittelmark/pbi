\name{pbi_etaSquared} 
\alias{pbi$etaSquared} 
\alias{pbi_etaSquared} 
\title{Effect size eta-squared for an Anova or a linear model object} 
\description{ 
    Cohen's rule of thumb for interpretation is: 0.01-0.09 small, 0.09-0.25 medium and above 0.25 large effect. You can convert eta-squared to a Pearson \emph{r} by using the square root of eta-squared. 
} 
\usage{pbi_etaSquared(x,y=NULL)} 
\arguments{ 
  \item{x}{ 
    vector with numerical values or a linear model or an aov object. 
  } 
  \item{y}{ 
    either f factor or NULL if x is given as model. 
  } 
} 
\value{return numerical value for eta-squares for all given variables in the model x or the value for the variable given in y if x is a numerical variable.} 
\examples{ 
  data(iris) 
  etaSquared=pbi_etaSquared 
  etaSquared(iris$Sepal.Length,iris$Species) 
  etaSquared(lm(iris$Sepal.Length ~ iris$Species)) 
  etaSquared(aov(iris$Sepal.Length ~ iris$Species)) 
  etaSquared(aov(Sepal.Length ~ Species+Sepal.Width,data=iris)) 
} 
\seealso{ \code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}},  
\code{\link[pbi:pbi_cohensD]{pbi_cohensD}}, \code{\link[pbi:pbi_epsilonSquared]{pbi_epsilonSquared}}} 
 
