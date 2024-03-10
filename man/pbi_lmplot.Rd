\name{pbi_lmplot} 
\alias{pbi$lmplot} 
\alias{pbi_lmplot} 
\title{ Plot a xy-plot with linear model fits and the confidence lines. } 
\description{ 
    The plot visualizes the linear fit description between numerical variables, 
  together with the confidence limits for the predictions as well as for 
  the linear model. The code is based on a tutorial by Conrad Halling (2006).  
  [https://web.archive.org/web/20180415155316/http://sphaerula.com/legacy/R/linearRegression.html](https://web.archive.org/web/20180415155316/http://sphaerula.com/legacy/R/linearRegression.html) 
} 
\usage{ pbi_lmplot(x,y,col="blue",pch=19,col.lm="red",grid=TRUE,...) } 
\arguments{ 
  \item{ x }{ 
    vector with numerical variables 
  } 
  \item{ y }{ 
    vector with numerical variables 
  } 
  \item{ col }{ 
    color for points and confidence lines for predictions, default: blue 
  } 
  \item{ pch }{ 
    plotting character, default: 19 
  } 
  \item{col.lm}{ 
    color for the linear model and the confidence lines, default: red 
  } 
  \item{ grid }{ 
    should a plot grid be drawn, default: TRUE 
  } 
  \item{\ldots}{other arguments delegated to the standard plot function} 
} 
\examples{ 
c20.22=c( 
17.9, 18.3, 18.3, 18.4, 18.4, 20.2, 20.3, 21.8, 21.9,  
22.1, 23.1, 24.2, 24.4) 
ins.sens=c( 
250, 220, 145, 115, 230, 200, 330, 400, 370, 260, 270,  
 530, 375) 
 
pbi_lmplot(x=c20.22,y=ins.sens, 
   xlab='\%C20-22 Fatty Acids',ylim=c(0,600), 
   xlim=c(17,25),main='best fit', 
   ylab='Insuline Sensitivity Index (mg/m^2/min)') 
 
legend('bottomright',c('best fit','fit 95\% CI', 
'prediction 95\% CI'),lty=c(1,2,2), 
col=c('red','red','blue')) 
} 
\seealso{  See also: [pbi](#pbi), [pbi_modelQuality](#modelQuality) } 
