\name{pbi_pcor.test} 
\alias{pbi$pcor.test} 
\alias{pbi_pcor.test} 
\title{Partial correlation test for two variables.} 
\description{ 
    Calculate partial correlation coefficient and  parametric  
  ("Pearson") or non-parametric ("Spearman")  
  test statistics for two variables corrected  
  for one or more other variables. 
} 
\usage{pbi_pcor.test(x,y,z,method='pearson')} 
\arguments{ 
  \item{x}{ 
    numeric vector, missing values are allowed 
  } 
  \item{y}{ 
    numeric vector, missing values are allowed 
  } 
  \item{z}{ 
    numeric vector, matrix or data frame,  missing values are allowed 
  } 
  \item{method}{ 
    character string indicating which partial correlation coefficient is to be computed, either "pearson" (default), or "spearman" 
  } 
} 
\value{return list with the following components:  
 
> - _estimate_ - gives the partial correlation coefficient between x and y given z 
  - _p.value_ - gives the p-value of the test 
  - _statistics_ - gives the value of the test statistics 
  - _n_ - gives the number of samples after deleting all the missing samples 
  - _gn_ - gives the number of given variables 
  - _method_ - gives the correlation method used} 
   
\examples{ 
  y.data = data.frame( 
   hl=c(7,15,19,15,21,22,57,15,20,18), 
   disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011), 
   deg=c(9,2,3,4,1,3,1,3,6,1), 
    BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00, 
         4.48e-03,2.10e-06,0.00e+00) 
  ) 
  # partial correlation between "hl" and "disp" given "deg" and "BC" 
  pbi_pcor.test(y.data$hl,y.data$disp,y.data[,c("deg","BC")]) 
} 
\seealso{  See also [pcor](#pcor) } 
