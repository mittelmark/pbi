\name{pbi_modelQuality} 
\alias{pbi$modelQuality} 
\alias{pbi_modelQuality} 
\title{Quality measures for regression and classification models} 
\description{ 
  The functions returns a few basic model quality measures  
  such as RMSE, MSE and MAE in case the data in x and y are numeric,  
  or Accuracy, Sensitivity, Specificity and Balanced Classification rate (BCR) 
  in case the data are factors. 
} 
\usage{ pbi_modelQuality(x,y) } 
\arguments{ 
  \item{ x }{ 
    numercial values or factor values representing the true values. 
  } 
  \item{ y }{ 
    numerical values or factor values representing the predicted values. 
  } 
} 
\value{return list with the components MAE, MSE, RMSE for numerical values or ACC, SENS, SPEC and BCR for factors (classes).} 
\examples{ 
library(MASS) 
data(swiss) 
# regression measures 
mod=lm(bwt~lwt,data=birthwt) 
summary(mod)$r.squared 
pred=predict(mod,newdata=birthwt) 
unlist(pbi_modelQuality(birthwt$bwt,pred)) 
mod=lm(bwt~lwt+age,data=birthwt) 
summary(mod)$r.squared 
unlist(pbi_modelQuality(birthwt$bwt,predict(mod,newdata=birthwt))) 
 
# classification measures 
library(rpart) 
# simple two class problem 
sex=as.factor(c(rep("M",50),rep("F",50))) 
height=c(rnorm(50,mean=180,sd=7),rnorm(50,mean=178,sd=6)) 
mod=rpart(sex~height) 
table(sex,pred=predict(mod,newdata=data.frame(height=height),type="class")) 
unlist(pbi_modelQuality(sex, 
       predict(mod,newdata=data.frame(height=height),type="class"))) 
# simple two class problem but unbalanced 
sex=as.factor(c(rep("M",30),rep("F",70))) # rare males 
height=c(rnorm(30,mean=180,sd=7),rnorm(70,mean=178,sd=6)) 
mod=rpart(sex~height) 
table(sex,pred=predict(mod,newdata=data.frame(height=height),type="class")) 
unlist(pbi_modelQuality(sex, 
       predict(mod,newdata=data.frame(height=height),type="class"))) 
# difficult three class problem 
table(iris$Species) 
mod=rpart(Species~.,data=iris) 
pred=predict(mod,newdata=iris,type="class") 
table(pred,iris$Species) 
# convert to two classes (setosa) 
table(c('s','v','v')[as.numeric(iris$Species)]) 
unlist(pbi_modelQuality( 
   as.factor(c('s','v','v')[as.numeric(iris$Species)]), 
   as.factor(c("s","v","v")[as.numeric(pred)]))) 
# convert to two classes (versicolor) 
unlist(pbi_modelQuality( 
   as.factor(c('x','v','x')[as.numeric(iris$Species)]), 
   as.factor(c("x","v","x")[as.numeric(pred)]))) 
# convert to two classes (virginica) 
unlist(pbi_modelQuality( 
   as.factor(c('x','x','v')[as.numeric(iris$Species)]), 
   as.factor(c("x","x","v")[as.numeric(pred)]))) 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}}, 
\code{\link[pbi:pbi_lmplot]{pbi_lmplot}} 
} 
