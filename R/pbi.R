## function documentation is directly embedded within
## the R code
## you can automatically extract the documentation using the bin/rman.R script
## by executing: using this like Rscript bin/rman.R R/add.R
## or from within the R-console:
## > setwd("package-folder")
## > source('bin/rman.R')
## > extractRd(list.files("R",pattern="*.R$",full.names=TRUE))

## package-documentation
#' \name{pbi-package}
#' \alias{pbi-package}
#' \title{pbi package - package template for R packages}
#' \description{The pbi package can be used as a template to create new packages from scratch.}
#' \details{Some more details:
#' The following list of objects and/or functions are available:
#' \describe{
#' \item{\link[pbi:add]{add(x,y)}}{an illustrative add function}
#' }
#' }
#' \examples{
#' library(pbi)
#' add(1,2)
#' }
#' \author{Detlef Groth <dgroth@uni-potsdam.de>}
#' \references{
#' \itemize{
#'    \item XYZ 
#' }
#' }
""

## Functions documentation, protect percent signs % with backslashes \%

#' \name{pbi_clusterSilhouette}
#' \alias{pbi$clusterSilhouette}
#' \alias{pbi_clusterSilhouette}
#' \title{Determine the cluster silhouette for given cluster ids and distance matrix.}
#' \description{
#'     This function can be used to determine the strength of a given clustering.
#'     A silhouette index of larger than 0.7 indicates a strong, 
#'     values between 0.5 and 0.7 a reasonable, between 0.25 and 0.5 a weak and below
#'     0.25 no reasonable structure.
#' }
#' \usage{pbi_clusterSilhouette(x,D)}
#' \arguments{
#'   \item{x}{cluster ids}
#'   \item{D}{a distance matrix object}
#' }
#' \value{return a list with the components avg.width (the silhouette index value for the complete) and clus.avg.widths the indices for the individual clusters.}
#' \examples{
#'   D=dist(scale(iris[,1:4]))
#'   hcl=hclust(D)
#'   pbi_clusterSilhouette(cutree(hcl,2),D)
#'   pbi_clusterSilhouette(cutree(hcl,3),D)$avg.width
#'   pbi_clusterSilhouette(cutree(hcl,4),D)$avg.width
#'   pbi_clusterSilhouette(cutree(hcl,5),D)$avg.width
#'   # three cluster seems to be the best
#' }

pbi_clusterSilhouette = function (x,D) {
  res=cluster::silhouette(x,D)
  sres=summary(res)
  return(list(avg.width=sres$avg.width,
              clus.avg.widths=sres$clus.avg.widths))
}

#' \name{pbi_clusterSimIndex}
#' \alias{pbi$clusterSimIndex}
#' \alias{pbi_clusterSimIndex}
#' \title{Compute similarity indices of two clusterings.}
#' \description{
#'   This function computes the similarity indices between two clusterings, 
#'   such as the Rand, Jaccard and Cohen's Kappa index.
#' }
#' \usage{pbi_clusterSimIndex(v1, v2)}
#' \arguments{
#'   \item{v1}{
#'     vector with cluster ids for first clustering
#'   }
#'   \item{v2}{
#'     vector with cluster ids for second clustering
#'   }
#' }
#' \value{return list with the components ...}
#' \examples{
#'   set.seed(123)
#'   hcl1=hclust(dist(iris[,1:4]))
#'   hcl2=hclust(dist(scale(iris[,1:4])))
#'   round(unlist(pbi_clusterSimIndex(cutree(hcl1,3),
#'     cutree(hcl2,3))),2)
#'   round(unlist(pbi_clusterSimIndex(cutree(hcl2,3),
#'     cutree(hcl1,3))),2) # symmetric
#'   round(unlist(pbi_clusterSimIndex(cutree(hcl1,3),
#'     sample(cutree(hcl2,3)))),2) # random values I
#'   round(unlist(pbi_clusterSimIndex(sample(cutree(hcl1,3)),
#'     sample(cutree(hcl2,3)))),2) # random values II
#' }

pbi_clusterSimIndex = function(v1, v2){
  ARI = function (cl1,cl2) {
    N=length(cl1)
    CM=table(cl1,cl2)
    a=apply(CM,1,sum)
    b=apply(CM,2,sum)
    SumNji=0
    SumCombai=0;
    SumCombbj=0;
    for (i in 1:max(cl1)) {
      for (j in 1:max(cl2)) {
        if (CM[i,j]>1) {
          SumNji=SumNji + choose(CM[i,j],2)
        }
      }
    }
    for (i in 1:max(cl1)) {
      if (a[i]>1) {
        SumCombai=  SumCombai+choose(a[i],2)
      }
    }
    for (j in 1:max(cl2)) {
      if ( b[j]>1) {
        SumCombbj=SumCombbj+choose(b[j],2)
      }
    }
    nCh2=choose(N,2);
    temp=(SumCombai*SumCombbj)/nCh2;
    AR =as.numeric((SumNji-temp)/(0.5*(SumCombai+SumCombbj)-temp))
    return(AR)
  }
  
  SS <- 0
  DD <- 0
  SD <- 0
  DS <- 0
  for (i in 1:(length(v1)-1)){
    for (j in (i+1):length(v1)){
      if (v1[i]==v1[j]){
        if (v2[i]==v2[j]){
          SS <- SS+1
        } else {
          SD <- SD+1
        }
      } else {
        if (v2[i]==v2[j]){
          DS <- DS+1
        } else {
          DD <- DD+1
        }
      }
    }
  }
  all = SS+DD+SD+DS
  po = (SS+DD)/all
  pe = (((SS+SD)/all)*((SS+DS)/all)) + (((DS+DD)/all)*((SD+DD)/all))
  Kappa = 1 - ((1-po) / (1-pe)) # same value as ARI!
  m <- min(max(v1), max(v2))
  tab=table(v1,v2)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  p0 <- sum(diag(tab[1:m, 1:m]))/n
  pc <- sum((ni[1:m]/n) * (nj[1:m]/n))    
  Kappa=(p0 - pc)/(1 - pc)
  # Mismatch!
  RI <- list(SS=SS, DD=DD, SD=SD, DS=DS, 
             Rand=(SS+DD)/(SS+DD+SD+DS),AdjRand=ARI(v1,v2),
             Jaccard=SS/(SS+SD+DS))
  return (RI)
}

#' \name{pbi_cohensD}
#' \alias{pbi$cohensD}
#' \alias{pbi_cohensD}
#' \title{ Effect size for the difference between two means.}
#' \description{
#'   The function *pbi_cohensD* calculates the effect size for the difference between two means.
#'   Due to Cohen's rule of thumb values of 0.2 to 0.5 are considered to stand 
#'   for small effects, values from 0.5 to 0.8 represent medium effects and values above 0.8 represent large effects.
#' }
#' \usage{pbi_cohensD(num,cat,paired=FALSE)}
#' \arguments{
#'   \item{num}{
#'     vector with numerical values
#'   }
#'   \item{cat}{
#'     vector with two grouping variables, having the same length as num
#'   }
#'   \item{paired}{
#'     are the data paired, default: FALSE
#'   }
#' }
#' \value{return Cohen's d value}
#' \examples{
#'   set.seed(125)
#'   data(sleep)
#'   with(sleep,pbi_cohensD(extra,group))
#'   x1=rnorm(100,mean=20,sd=1)
#'   x2=rnorm(100,mean=22,sd=1)
#'   g1=rep('A',100)
#'   g2=rep('B',100)
#'   # difference should be around 2SD
#'   pbi_cohensD(c(x1,x2),as.factor(c(g1,g2)))
#'   # biserial correlation coefficient as alternative
#'   # value is as well large
#'   cor(c(x1,x2),as.numeric(as.factor(c(g1,g2))))
#' }
#' \seealso{  See also [pbi](#home), [pbi_cohensW](#cohensW) }

pbi_cohensD <- function (num, cat,paired=FALSE) {
  if (paired) {
    tt=t.test(num ~ cat,paired=paired)
    return(tt$statistic[[1]]/sqrt(length(num/2)))
  }   
  tt.agg=aggregate(num,by=list(cat),
                   mean,na.rm=TRUE)
  pooledSD <- function(x, y) {
    x=x[!is.na(x)]
    y=y[!is.na(y)]
    sq.devs <- (c(x - mean(x), y - mean(y)))^2
    n <- length(sq.devs)
    return(sqrt(sum(sq.devs)/(n - 2)))
  }
  d=abs(tt.agg$x[1]-tt.agg$x[2])/pooledSD(
    num[cat==levels(cat)[1]],
    num[cat==levels(cat)[2]])
  return(d)
}

#' \name{pbi_cohensH}
#' \alias{pbi$cohensH}
#' \alias{pbi_cohensH}
#' \title{Effect size for a 2x2 contingency table.}
#' \description{
#'   The function *pbi_cohensH* calculates the effect size for 2x2 contingency tables. 
#'   Due to Cohen's rule of thumb values of 0.2 to 0.5 are considered to stand 
#'   for small effects, values from 0.5 to 0.8 represent medium effects and values above 0.8 represent large effects. 
#' }
#' \usage{pbi_cohensH(tab)}
#' \arguments{
#'   \item{tab}{
#'     a 2x2 contingency table
#'   }
#' }
#' \value{return Cohen's h value}
#' \examples{
#' # data from New Eng. J. Med. 329:297-303, 1993
#' azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2))
#' rownames(azt)=c("AZT","Placebo")
#' colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress")
#' pbi_cohensH(azt)
#' }
#' \seealso{  See also [pbi](#home), [pbi_cohensW](#cohensW) }

pbi_cohensH <- function (tab) {
  pt=prop.test(tab)
  h=2*abs(asin(sqrt(pt$estimate[1]))-
            asin(sqrt(pt$estimate[2])))
  return(h[[1]])
}

#' \name{pbi_cohensW}
#' \alias{pbi$cohensW}
#' \alias{pbi_cohensW}
#' \title{Effect size for 2x2 and larger contingency tables as well as for single variables.}
#' \description{
#'   The function *pbi_cohensW* calculates the effect size for contingency tables. 
#'   Due to Cohen's rule of thumb values of 0.1 to 0.3 are considered to stand 
#'   for small effects, values from 0.3 to 0.5 represent medium effects and values 
#'   above 0.5 represent large effects.
#' }
#' \usage{pbi_cohensW(x,p=NULL)}
#' \arguments{
#'   \item{x}{
#'     contingency table with counts, usually created using the table 
#'     command for two variables  or vector with counts with observations
#'     for a single variable, if x is a vector p must be given
#'   }
#'   \item{p}{
#'     expected proportions if x is a vector, in case *x* has length of two,
#'     a single value can be given, otherwise p must have same length as x, default: NULL
#'   }
#' }
#' \examples{
#' data(Titanic)
#' Titanic[1,1,,]
#' pbi_cohensW(Titanic[1,1,,])
#' # Data from New Eng. J. Med. 329:297-303, 1993
#' azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2))
#' rownames(azt)=c("AZT","Placebo")
#' colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress")
#' pbi_cohensW(azt)
#' # number of boys (25) and girls (15) in a hospital which deviates 
#' # from 50/50 or 51/49 ratios
#' pbi_cohensW(c(25,15),p=0.5)
#' pbi_cohensW(c(25,15),p=c(0.51,0.49))
#' # most extrem case 40 boys and 0 girls
#' pbi_cohensW(c(40,0),p=0.5)
#' pbi_cohensW(c(40,0),p=c(0.51,0.49)) # max value here around 2*0.49
#' }
#' \seealso{  [pbi](#home), [pbi_cohensH](#cohensH) }

pbi_cohensW <- function (x,p=NULL) {
  if (is.table(x) | is.matrix(x)) {
    tab=x
    pe=prop.table(chisq.test(tab)$expected)
    po=prop.table(tab)
    w=sqrt(sum(((po-pe)^2)/pe))
    return(w[[1]])
  } else if (is.null(p)) {
    stop('Error: If x is a vector, p must be given!')
  } else {
    if (length(x) == 2 & length(p) == 1) {
      p=c(p,1-p)
      po=prop.table(x)
      pe=p
    } else if  (length(x) == length(p)) {
      po=prop.table(x)
      pe=p
    } else {
      stop('Error: for more than 2 categories the
                 given proportion vector p must have the
                 same length as the given count vector x')
    }
    w = sqrt(sum(((po-pe)^2)/pe))
    return(w)
  }
}

#' \name{pbi_cv}
#' \alias{pbi$cv}
#' \alias{pbi_cv}
#' \title{Calculate the coefficient of variation.}
#' \description{
#'   The function *pbi_cv* calculates the coefficient of variation, which is
#'   a unitless measure of data scatter. It is calculated as: _(100*sd(x))/mean(x)_.
#'   All data values has to be non-negative.
#' }
#' \usage{pbi_cv(x,na.rm=FALSE)}
#' \arguments{
#'   \item{x}{
#'     vector with positive numerical values.
#'   }
#'   \item{na.rm}{
#'     should NA's be removed, default: FALSE
#'   }
#' }
#' \value{return numerical value for the coefficient of variation.}
#' \examples{
#'   pbi_cv(rnorm(20,mean=100,sd=4))
#'   pbi_cv(c(1,2,3,4))
#'   pbi_cv(c(1,2,3,4,NA))
#'   pbi_cv(c(1,2,3,4,NA),na.rm=TRUE)
#' }

pbi_cv = function (x,na.rm=FALSE) {
  cv=100*sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm)
  return(cv)
}

#' \name{pbi_dassoc}
#' \alias{pbi$dassoc}
#' \alias{pbi_dassoc}
#' \title{Create assocplots with residual coloring.}
#' \description{
#'    This function updates the standard assocplot function from the graphics package 
#'    with the ability to display residual colors. In blue and red are shown groups with 
#'    residuals above +4 or below -4 in light colors are shown residuals between 2 and 4 for positive and -4 and -2 for negative residuals.
#' }
#' \usage{pbi_dassoc(...,shade=TRUE)}
#' \arguments{
#'   \item{...}{
#'     arguments delegated to the standard assocplot function.
#'   }
#'   \item{shade}{
#'     should the residuals been shown, default: TRUE
#'   }
#' }
#' \examples{
#'   x = margin.table(HairEyeColor, c(1, 2))
#'   pbi_dassoc(x)
#' }

pbi_dassoc <- function (...,shade=TRUE) {
  # https://stackoverflow.com/questions/38732663/how-to-insert-expression-into-the-body-of-a-function-in-r
  funins <- function(f, expr = expression(x<-2*x), after=1) {
    body(f)<-as.call(append(as.list(body(f)), expr, after=after))
    f
  }
  
  expr=expression({
    # DG: changed for shade
    residuals=chisq.test(x)$residuals
    cols=c('#CF3761','#E18E9E','#E0E0E0','#96A2DF','#4267E0')
    resis=c()
    # R plots from top right to lower left ...
    # we rearrange the colors
    for (c in ncol(x):1) {
      resis=c(resis,residuals[,c])
    }
    acols=cols[cut(resis,
                   breaks=c(-Inf,-4,-2,2,4,Inf),
                   labels=c(1,2,3,4,5))]
    rect(z[, 1] - e/2, z[, 2], z[, 1] + e/2, z[, 2] + d, col = acols)
  })
  if (shade) {
    g <- funins(graphics::assocplot, expr,after=length(body(graphics::assocplot)))
    g(...)
  } else {
    graphics::assocplot(...)
  }
  
}

#' \name{pbi_dcorr}
#' \alias{pbi$dcorr}
#' \alias{pbi_dcorr}
#' \title{Calculate pairwise correlations for a given data frame or matrix including their p-values.}
#' \description{
#'   The function is an extension to the standard _stats::cor_ function, it calculates as well
#'   the p-values for the pairwise aoosciations and returns them in a matrix as well.
#' }
#' \usage{pbi_dcorr(data,method='pearson',use='pairwise.complete.ob')}
#' \arguments{
#'   \item{data}{
#'     matrix or data frame where the variables are in the columns, NA's are allowed.
#'   }
#'   \item{method}{
#'     type of correlation to be determined, either 'pearson', 'spearman' or 'kendall', default: 'pearson'
#'   }
#'   \item{use}{
#'     how to deal with NA's, default: 'pairwise.complete.obs'
#'   }
#' }
#' \value{return list with the following components:
#' 
#' > - _estimate_ - matrix with correlation values
#'   - _p.value_ - matrix with p-values
#'   - _method_ - character string with the used correlation method
#'   }
#' \examples{
#'   data(swiss)
#'   lapply(pbi_dcorr(swiss)[1:2],round,2)
#' }
#' \seealso{  See also [pbi](#home), [pbi_dcorrplot](#dcorrplot) }

pbi_dcorr <- function (data,method='pearson',use='pairwise.complete.ob') {
  mt=matrix(0,nrow=ncol(data),ncol=ncol(data))
  colnames(mt)=rownames(mt)=colnames(data)
  mt.pval=mt
  diag(mt)=1
  for (i in 1:(ncol(data)-1)) {
    for (j in i:ncol(data)) {
      rt=cor.test(data[,i],data[,j],
                  method=method,use=use)
      mt[i,j]=mt[j,i]=rt$estimate
      mt.pval[i,j]=mt.pval[j,i]=rt$p.value
    }
  }
  return(list(estimate=mt,p.value=mt.pval,method=method))
}

#' \name{pbi_dcorrplot}
#' \alias{pbi$dcorrplot}
#' \alias{pbi_dcorrplot}
#' \title{Visualize a correlation matrix.}
#' \description{
#'     Visualize a correlation matrix.
#' }
#' \usage{pbi_dcorrplot(mt,text.lower=TRUE, text.upper=FALSE,pch=19,p.mat=NULL,alpha=0.05,cex.sym=5,cex.r=1,cex.lab=1.4,...)}
#' \arguments{
#'   \item{mt}{
#'     matrix with pairwise correlations
#'   }
#'   \item{text.lower}{
#'     should in the lower diagonal the correlation coefficient be shown, default: TRUE
#'   }
#'   \item{text.upper}{
#'     should in the upper diagonal the correlation coefficient be shown, default: FALSE
#'   }
#'   \item{pch}{
#'     the plotting symbol for the correlations, default: 19
#'   }
#'   \item{p.mat}{
#'     matrix with p-values to strike out insignificant p-values, default: NULL (not used)
#'   }
#'   \item{alpha}{
#'     significance threshold for _p.mat_, default: 0.05
#'   }
#'   \item{cex.sym}{
#'     character expansion for the correlation symbols, default: 5
#'   }
#'   \item{cex.r}{
#'     character expansion for the r-values if _text.lower_ or _text.upper_ are set to TRUE, default: 1
#'   }
#'   \item{cex.lab}{
#'     character expansion for the variable text labels, default: 1.4
#'   }
#'   \item{...}{
#'     arguments delegated to the plot function
#'   }
#' }
#' \examples{
#'   data(swiss)
#'   sw=swiss
#'   colnames(sw)=abbreviate(colnames(swiss),6)
#'   cr=pbi_dcorr(sw,method='spearman')
#'   pbi_dcorrplot(cr$estimate,cex.sym=8,text.lower=TRUE,
#'      cex.r=1.5,p.mat=cr$p.value)
#' }
#' \seealso{  See also [pbi](#home), [pbi_dcorr](#dcorr) }

pbi_dcorrplot <- function (mt,text.lower=TRUE, text.upper=FALSE,
                             pch=19,p.mat=NULL,alpha=0.05,
                             cex.sym=5,cex.r=1,cex.lab=1.4,...) {
  if (!is.matrix(p.mat)) {
    p.mat=mt
    p.mat[]=0
  }
  yend=nrow(mt)+1
  xend=ncol(mt)+1
  plot(1,type="n",xlab="",ylab="",axes=FALSE,
       xlim=c(-0.5,xend),ylim=c(nrow(mt)+0.35,0),...)
  text(1:(ncol(mt)),0.25,colnames(mt),cex=cex.lab)
  text(-0.5,1:nrow(mt),rownames(mt),cex=cex.lab,pos=4)
  cols=paste("#DD3333",rev(c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD")),sep="")
  cols=c(cols,paste("#3333DD",c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD"),sep=""))
  breaks=seq(-1,1,by=0.1)                  
  sym=identical(rownames(mt),colnames(mt))
  for (i in 1:nrow(mt)) {
    for (j in 1:nrow(mt)) {
      if (sym & i == j) {
        next
      }   
      coli=cut(mt[i,j],breaks=breaks,labels=1:20)
      if (i == j & !sym & text.lower) {
        text(i,j,round(mt[i,j],2),cex=cex.r)
        if (p.mat[i,j]>alpha) {
          text(i,j,"x",cex=cex.r*2)
        }
        
      } else if (i < j & text.lower) {
        text(i,j,round(mt[i,j],2),cex=cex.r)
        if (p.mat[i,j]>alpha) {
          text(i,j,"x",cex=cex.r*2)
        }
        
      } else if (i > j & text.upper) {
        text(i,j,round(mt[i,j],2),cex=cex.r)
        if (p.mat[i,j]>alpha) {
          text(i,j,"x",cex=cex.r*2)
        }
      } else {
        points(i,j,pch=pch,cex=cex.sym,col=cols[coli])
        if (p.mat[i,j]>alpha) {
          text(i,j,"x",cex=cex.sym*0.3)
        }
      }
    }
  }
}

#' \name{pbi_df2md}
#' \alias{pbi$df2md}
#' \alias{pbi_df2md}
#' \title{Convert a data frame or a matrix into a Markdown table.}
#' \description{
#'   This function can be used within Rmarkdown documents to display easily
#'   a simple Markdown table. For more advance use  cases you should other commands
#'   such as kable from the knitr package.
#' }
#' \usage{pbi_df2md(x,caption='',rownames=TRUE)}
#' \arguments{
#'   \item{x}{
#'     matrix or data frame 
#'   }
#'   \item{rownames}{
#'     should  the rownames be displayed, default: TRUE
#'   }
#'   \item{caption}{
#'     the caption for the table, it is just displayed below of the table.
#'   }
#' }
#' \value{return prints to stdout.}
#' \examples{
#'   data(swiss)
#'   pbi_df2md(head(swiss))
#' }

pbi_df2md <- function(x,caption='',rownames=TRUE) {
  df=x
  cn <- colnames(df)
  if (is.null(cn[1])) {
    cn=as.character(1:ncol(df))
  }
  rn <- rownames(df)
  if (is.null(rn[1])) {
    rn=as.character(1:nrow(df))
  }
  if (rownames) {
    headr <- paste0(c("","", cn),  sep = "|", collapse='')
    sepr <- paste0(c('|', rep(paste0(c(rep('-',3), "|"), 
                                     collapse=''),length(cn)+1)), collapse ='')
  } else {
    headr <- paste0(c("", cn),  sep = "|", collapse='')
    sepr <- paste0(c('|', rep(paste0(c(rep('-',3), "|"), 
                                     collapse=''),length(cn))), collapse ='')
    
  }
  st <- "|"
  for (i in 1:nrow(df)){
    if (rownames) {
      st <- paste0(st, "**",as.character(rn[i]), "**|", collapse='')
    }
    for(j in 1:ncol(df)){
      if (j%%ncol(df) == 0) {
        st <- paste0(st, as.character(df[i,j]), "|", 
                     "\n", "" , "|", collapse = '')
      } else {
        st <- paste0(st, as.character(df[i,j]), "|", 
                     collapse = '')
      }
    }
  }
  fin <- paste0(c(headr, sepr, substr(st,1,nchar(st)-1)), collapse="\n")
  if (caption!='') {
    fin=paste0(fin,'\n',caption,'\n')
  }
  cat(fin)
}

#' \name{pbi_dist}
#' \alias{pbi$dist}
#' \alias{pbi_dist}
#' \title{Distance function which as well supports binary and correlation distances in addition to standard distances.}
#' \description{
#'   This function is an extension to stats::dist function as  it supports as well correlation distance and 
#'   binary distance measures such as Jaccard coefficient and Matching coefficient. 
#'   The correlation distance is implemented as `D(i,j) = 1-((r(i,j)+1)/2)` 
#'   so negative correlations have low similarities.
#' }
#' \usage{pbi_dist(x,method="euclidean",...)}
#' \arguments{
#'   \item{x}{
#'     data frame or matrix with numerical data, 
#'     in case of binary data as well boolean and two level nominal data could be supplied.
#'   }
#'   \item{method}{
#'     the distance measure to be used, one of the measures 
#'     for `stats::dist` such as "euclidean" or "correlation", 
#'     alias for  "pearson" or "spearman","kendall", 
#'     for binary data "jc" (Jaccard) and "mc" (Matching coeffient) are supported.
#'   }
#'   \item{...}{
#'     remaining arguments are forwarded to stats::dist function in case the method is 
#'     handled by this default method. 
#'   }
#' }
#' \value{return distance matrix}
#' \examples{
#'   round(cor(iris[,1:4]),2)
#'   round(as.matrix(pbi_dist(t(iris[,1:4]),method="pearson")),2)
#'   biris=iris[,1:4]
#'   biris=apply(biris,2,function(x) { return(x>median(x)) })
#'   head(biris,3)
#'   head(apply(biris,2,as.numeric),3)
#'   biris=apply(biris,2,as.numeric)
#'   summary(biris)
#'   round(as.matrix(pbi_dist(t(biris),method="mc")),2)
#'   d.can=as.matrix(pbi_dist(t(scale(iris[,1:4])),method="can"))
#'   d.can
#'   d.can=d.can/max(d.can)
#'   round(d.can,2)
#' }

pbi_dist = function (x,method="euclidean",...) {
  bdist <- function (mt,method='mc') {
    n = ncol(mt)
    res = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        freq = mt[,i]+mt[,j]
        d = length(which(freq == 0))
        a = length(which(freq == 2))
        bc = length(which(freq == 1))
        if (method=='mc') {
          res[i,j] = res[j,i] = (a + d)/(a + d + bc)
        } else if (method == 'jc') {
          res[i,j] = res[j,i] = a/(a + bc) 
        } else {
          stop('unknown method, known methods are 
                         jc - jaccard and mc matching coefficient')
        }
      }
    }   
    rownames(res)=colnames(mt)
    colnames(res)=colnames(mt)
    res=1-res
    return(as.dist(res))
  }
  if (method %in% c("correlation","cor","pearson","spearman","kendall")) {
    if (method %in% c("correlation","cor")) {
      method = "pearson"
    }
    D=as.dist(1-((cor(t(x),method=method,use="pairwise.complete.obs")+1)/2))
    return(D)
  } else if (method %in% c("mc","jc")) {
    return(bdist(t(x),method=method))
  } else {
    return(stats::dist(x,method=method,...))
  }
}

#' \name{pbi_domainPlot}
#' \alias{pbi$domainPlot}
#' \alias{pbi_domainPlot}
#' \title{Use the mydomain website of prosite to create a protein domain plot.}
#' \description{
#'   This function can be used to draw a plot representing protein domains into a standard R plot.
#'   The data are send to the website [https://prosite.expasy.org/cgi-bin/prosite/mydomains/](https://prosite.expasy.org/cgi-bin/prosite/mydomains/), a URL is created and the image from this url is downloaded to the local file system. The filename is a CRC32 digest of the URL so allowing to cache as well the downloaded results.
#' }
#' \usage{pbi_domainPlot(domains,ranges,sites,length=1000,hscale=1,cache=TRUE,plot=TRUE,cex=0.8,...)}
#' \arguments{
#'   \item{domains}{
#'     list where names are the domain names and the values are four integers: start, end, shape (1:6), color (1:4)
#'   }
#'   \item{ranges}{
#'     list with vectors of three integers: start, end, type (0:1)
#'   }
#'   \item{sites}{
#'     list with vectors of two integers: position, type (0:1)
#'   }
#'   \item{length}{
#'     length of the protein, default: 1000
#'   }
#'   \item{hscale}{
#'     horizontal scaling factor (between 0.1 and 2.0), default: 1.0 
#'   }
#'   \item{cache}{
#'     should the image locally cached, default: TRUE
#'   }
#'   \item{plot}{
#'     should the image plotted on a R plot devices, default: TRUE
#'   }
#'   \item{cex}{
#'     character expansion for the domain text and sites if plot is TRUE, default: 0.8
#'   }
#'   \item{...}{
#'     any argument forwarded to the plot function if plot is TRUE
#'   }
#' }
#' \value{return a URL or a local filename in case of cache=TRUE, if plot is TRUE NULL is returned invisible}
#' \examples{
#' url=pbi_domainPlot(
#'     domains=list(MYDOM1=c(100,200,2,1),MYDOM2=c(300,500,3,2)),
#'     ranges=list(a=c(190,500,1)),
#'     sites=list(a=c(150,1)),hscale=2.0,plot=FALSE)
#' print(url)
#' }

pbi_domainPlot <- function (domains,ranges,sites,length=1000,hscale=1,cache=TRUE,plot=TRUE,cex=0.8,...) {
  shape <- function (x,y,width=1,height=0.3,type="circle",arrow=TRUE,dir="left") {
    center = function (poly) {
      poly$x=poly$x-mean(poly$x)
      poly$x=poly$x/diff(range(poly$x))
      poly$y=poly$y-mean(poly$y)
      poly$y=poly$y/diff(range(poly$y))
      return(poly)
    }
    circle = function (x,y, radius=1,length=100) {
      theta = seq(0, 2 * pi, length = 100) 
      return(list(x=radius*cos(theta)+x,
                  y=radius*sin(theta)+y))
    }
    if (substr(type,1,4) == "circ") {
      poly=circle(0,0,radius=1)
    } else if (type %in% c("rect","rectangle")) {
      poly=list(x=c(-0.5,0.5,0.5,-0.5),
                y=c(-0.5,-0.5,0.5,0.5))
    } else if (type == "diamond") {
      poly=list(x=c(-0.5,  0,   0.5,0  ),
                y=c( 0  , -0.5, 0.0,0.5))
    } else if (type == "hexagon") {
      poly=list(x=c(-0.5,-0.3,+0.3,+0.5,+0.3,-0.3),
                y=c(+0.0,-0.5,-0.5,+0.0,+0.5,+0.5))
      if (arrow) {
        poly$x[1]=-0.05
      }
      if (dir == "right") {
        poly$x=poly$x*-1
      }
      if (dir == "top") {
        x=poly$x
        poly$x=poly$y
        poly$y=x
      }
    } else if (type == "octagon") {
      poly=list(x=c(-0.5,-0.5,-0.3,0.3,0.5,0.5,0.3,-0.3),
                y=c(0.25,-0.25,-0.5,-0.5,-0.25,+0.25,0.5,0.5))
    } else if (type == "pentagon") {
      poly=list(x=c(-0.5,-0.25,0.5,0.5,-0.25),
                y=c(0,0.5,0.5,-0.5,-0.5))
      if (arrow) {
        poly$y[3]=0.25
        poly$y[4]=-0.25
      }
      if (dir == "right") {
        poly$x=poly$x*-1
      }
    } else {
      stop("Error: Unkown type '",type,"'!",sep="")
    }
    poly=center(poly)
    poly$x=poly$x*width+x
    poly$y=poly$y*height+y
    return(poly)
  }
  if (plot) {
    plot(0,type="n",xlim=c(0,length),ylim=c(-0.3,0.5),axes=FALSE,xlab="",ylab="",...)
    lines(x=c(0,length),y=c(0,0),lwd=2)
    cols=c("orange","#99ff99","#3399ff","#aaaaaa",rainbow(20))
    shapes=list(rect=shape(x=0,y=0,type="rect"),
                circle=shape(x=0,y=0,type="circle"),
                pentag1=shape(x=0,y=0,type="pentagon",dir='right',arrow=TRUE),
                pentag2=shape(x=0,y=0,type="pentagon",dir='left',arrow=TRUE),
                hexagon1=shape(x=0,y=0,type="pentagon",dir='left'),
                hexagon2=shape(x=0,y=0,type="pentagon",dir='top'))
    for (range in ranges) {
      lines(x=c(range[1],range[2]),y=c(0.3,0.3),lwd=2)
      if (range[3] == 1) {
        lines(x=c(range[1],range[1]),y=c(0.3,0),lwd=2)
        lines(x=c(range[2],range[2]),y=c(0.3,0),lwd=2)
      }
    }
    for (site in sites) {
      col="grey30"
      if (site[2] == 1) {
        col="red"
      }
      lines(x=c(site[1],site[1]),y=c(0.2,0),lwd=2,col=col)
      points(site[1],0.2,pch=18,cex=cex*2,col=col)
    }
    for (domain in names(domains)) {
      vec=domains[[domain]]
      width=vec[2]-vec[1]
      x=(vec[2]+vec[1])/2
      shape=shapes[[vec[3]]]
      shape$x=shape$x*width
      shape$x=shape$x+x
      polygon(shape,col=cols[vec[4]])
      text(x,0,domain,cex=cex)
    }
  } else {
    url="https://prosite.expasy.org/cgi-bin/prosite/PSImage.cgi?"
    i=1
    for (domain in names(domains)) {
      vec=domains[[domain]]
      if (i>1) {
        url=paste(url,"+",sep="")
      } else {
        url=paste(url,"hit=",sep="")
      }
      url=sprintf("%s%i,%i,%i_%i,%s",url,vec[1],vec[2],vec[3],vec[4],domain)
      i=i+1
    }
    i=1
    for (range in ranges) {
      if (i>1) {
        url=paste(url,"+",sep="")
      } else {
        url=paste(url,"&range=",sep="")
      }
      i=i+1
      vec=range
      url=sprintf("%s%i,%i,%i",url,vec[1],vec[2],vec[3])
    }
    i=1
    for (site in sites) {
      if (i>1) {
        url=paste(url,"+",sep="")
      } else {
        url=paste(url,"&site=",sep="")
      }
      i=i+1
      vec=site
      url=sprintf("%s%i,%i",url,vec[1],vec[2])
    }
    url=paste(url,"&len=",length,"&scale=",hscale,sep="")
    if (!file.exists("img")) {
      dir.create("img")
    }
    if (cache) {
      if (!requireNamespace("digest",quietly=TRUE)) {
        stop("Error: caching requires library digest!")
      }
      dig=digest::digest(url,algo="crc32")
      file=file.path("img",paste(dig,".png",sep=""))
      if (!file.exists(file)) {
        download.file(url,file)
      }
      return(file)
    }
    return(url)
  }
}

#' \name{pbi_dpairs}
#' \alias{pbi$dpairs}
#' \alias{pbi_dpairs}
#' \title{Improved pairs plot considering the data types.}
#' \description{
#'   The function _dpairs_ provides an improved pairs plot which accounts
#'   for the data type of the actual variables. It will plot in the 
#'   lower diagonal xy-plots, box-plots or assoc-plots depending on the 
#'   two data types. In the upper diagonal effect sizes and stars for the p-values
#'   for the tests (anova, t.test, chisq.test or cor.test) will be shown. In the diagonal 
#'   the data distribution will be outlined. This plot is usually an useful visualization for 3-8 variables.
#' }
#' \usage{pbi_dpairs(data,col.box='grey80',col.xy="grey60",cex.diag=2.5,order=TRUE,pch=19)}
#' \arguments{
#'   \item{data}{
#'     data frame with columns of class factor, numeric or integer.
#'   }
#'   \item{col.box}{
#'     colors for the boxplots, either a single value or a vector of colors for each level of a factor variable, default; 'grey80'
#'   }
#'   \item{col.xy}{
#'     colors for the xy-plots, either a single value of a vector which is as long as the number of data points, default: 'grey60'
#'   }
#'   \item{cex.diag}{
#'     character expansion for the diagonal texts
#'   }
#'   \item{order}{
#'     should the variables be ordered by data type and name, this is recommended as it orders the plots, starting with assocplots, then boxplots and finally xyplots, default: TRUE
#'   }
#'   \item{pch}{
#'     plotting character for xy-plots, default 19 (round circle).
#'   }
#' }
#' \examples{
#'   data(iris)
#'   par(omi = c(0.8, 0.4,0.4,0.4))
#'   pbi_dpairs(iris,col.box=2:4,col.xy=rep(c(2:4),each=50),
#'      cex.diag=1.6)
#'   pbi_dpairs.legend(levels(iris$Species),col=2:4)
#' 
#'   par(omi=c(0.5,0.5,0.8,0.2))
#'   library(MASS)
#'   btwt=birthwt; 
#'   for (col in c('low','race','smoke','ptl','ht','ui','ftv')) { 
#'      btwt[,col]=as.factor(btwt[,col]) 
#'   }
#'   pbi_dpairs(btwt[,2:8],cex.diag=1.6)
#'       mtext('Birth-Weight data',side=3,outer=TRUE,
#'       cex=1.5,line=1)
#'   if (require("palmerpenguins")) {
#'       data(penguins)
#'       peng=penguins
#'       colnames(peng)[3]='bill.len'
#'       colnames(peng)[4]='bill.dep'
#'       colnames(peng)[5]='flip.len'
#'       colnames(peng)[6]='mass'
#'       pbi_dpairs(peng,col.xy=pbi_pastel(3)[as.numeric(peng$species)])
#'       pbi_dpairs.legend(levels(peng$species),
#'           col=pbi_pastel(3)[as.numeric(as.factor(levels(peng$species)))])
#'   }
#' }

pbi_dpairs <- function (data,col.box='grey80',col.xy="grey60",cex.diag=2.5,
                          order=TRUE,pch=19) {
  oop=options()
  options(warn=-1)
  report.pval=pbi_report.pval
  if (any(class(data) %in% "tbl_df")) {
    data=as.data.frame(data)
  }
  if (order) {
    data=data[,sort(colnames(data))]
    res=c(); for (i in 1:ncol(data)) { res=c(res,class(data[,i])) }
    idx=order(res)
    data=data[,idx]
  }
  mai=rep(0.0,4)
  opar=par(mfrow=c(ncol(data),ncol(data)),mai=mai)
  cnames=colnames(data)
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      if (i == j) {
        plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
        text(0.5,0.5,cnames[i],cex=cex.diag)
        par(mai=rep(0,4))
        box(lty=3,col='grey70')
        if (class(data[,i]) == "factor") {
          rect(0.1,0.1,0.9,0.3,col="grey90")
          for (ci in cumsum(prop.table(table(data[,i])))) {
            x=0.1+0.8*ci
            lines(x=c(x,x),y=c(0.1,0.3))
          }
        }
        if (class(data[,i]) %in% c("numeric","integer")) {
          ds=density(data[,i],na.rm=TRUE)
          ds$x=ds$x-min(ds$x)
          ds$x=ds$x/max(ds$x)
          ds$y=(ds$y/max(ds$y)*0.3)
          polygon(ds,col='grey80')
        }
        par(mai=mai)
      } else if (i > j) {
        if (class(data[,i]) %in% c("numeric","integer") & class(data[,j]) %in% c("numeric","integer")) {
          plot(data[,i] ~ data[,j],xlab='',ylab='',axes=FALSE,pch=pch,
               col=col.xy)
          box(col='grey70')
          if (j+1 == i) {
            #axis(3)
            #axis(4)
          }
          if (j == 1) {
            axis(2)
          }
          if (i == ncol(data)) {
            ticks=axTicks(1) 
            axis(1,at=ticks[1:(length(ticks)-1)],labels=ticks[1:(length(ticks)-1)],col='grey70')
          }
        } else if (class(data[,i]) == "factor" & class(data[,j]) == "factor") {
          par(mai=rep(0.3,4))
          pbi_dassoc(t(table(data[,i],data[,j])))
          par(mai=rep(0,4))
          box(lty=3,col='grey70')
          par(mai=mai)
        } else if (class(data[,i]) %in% c("numeric","integer")) {
          boxplot(data[,i] ~ data[,j],col=col.box,axes=FALSE)
          if (j+1 == i) {
            #axis(3,at=1:length(levels(data[,j])),labels=levels(data[,j]))
            #axis(4)
          } 
          if (j == 1) {
            ticks=axTicks(2) 
            axis(2,at=ticks[1:(length(ticks)-1)],labels=ticks[1:(length(ticks)-1)],col='grey70')
          }
          if (i == ncol(data)) {
            axis(1,at=1:length(levels(data[,j])),labels=levels(data[,j]),col="grey70")
          }
          
          box(col="grey70")
        } else if (class(data[,j]) %in% c("numeric","integer")) {
          boxplot(data[,j] ~ data[,i],col=col.box,axes=FALSE)
          if (j == 1) {
            axis(2)
          }
          if (i == ncol(data)) {
            axis(1,at=1:length(levels(data[,j])),labels=levels(data[,j]))
          }
          box()
          
        } 
      } else {
        if (class(data[,i]) %in% c("numeric","integer") & class(data[,j]) %in% c("numeric","integer")) {
          r=cor.test(data[,i],data[,j])
          rs=cor.test(data[,i],data[,j],method='spearman')
          plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
          text(0.5,0.59,bquote("" ~ r[P] ~ .(sprintf(" = %.2f%s",r$estimate,report.pval(r$p.value,star=TRUE)))),cex=1.5)
          text(0.5,0.41,bquote("" ~ r[S] ~ .(sprintf(" = %.2f%s",rs$estimate,report.pval(rs$p.value,star=TRUE)))),cex=1.5)
        } else if (class(data[,i]) == "factor" & class(data[,j]) == "factor") {
          cw=pbi_cohensW(table(data[,i],data[,j]))
          chsq=chisq.test(table(data[,i],data[,j]))
          plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
          text(0.5,0.5,sprintf("Cohen's w =\n%.2f %s",cw,report.pval(chsq$p.value,star=TRUE)),cex=1.5)
          
        } else if (class(data[,i]) %in% c("numeric","integer")) {
          if (length(levels(data[,j]))==2) {
            tt=t.test(data[,i] ~ data[,j]) 
            cd=pbi_cohensD(data[,i],data[,j])
            plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
            text(0.5,0.5,sprintf("Cohen's d =\n%.2f %s",cd,report.pval(tt$p.value,star=TRUE)),cex=1.5)
          } else {
            raov=aov(data[,i] ~ data[,j]) 
            #recover()
            rs=pbi_etaSquared(raov)
            pval=report.pval(summary(raov)[[1]][1,5],star=TRUE)
            plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
            text(0.5,0.5,bquote(eta~2~sprintf(" = %.2f %s",rs,pval)),cex=1.5)
          }
        } else if (class(data[,j]) %in% c("numeric","integer")) {
          if (length(levels(data[,i]))==2) {
            tt=t.test(data[,j] ~ data[,i]) 
            cd=pbi_cohensD(data[,j],data[,i])
            plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
            text(0.5,0.5,sprintf("Cohen's d =\n%.2f %s",cd,report.pval(tt$p.value,star=TRUE)),cex=1.5)
          } else {
            raov=aov(data[,j] ~ data[,i]) 
            rs=pbi_etaSquared(raov)
            pval=report.pval(summary(raov)[[1]][1,5],star=TRUE)
            plot(1,type='n',xlab='',ylab='',axes=FALSE,xlim=c(0,1),ylim=c(0,1))
            val=sprintf("%.2f %s",rs,pval)
            text(0.5,0.5,bquote("" ~ eta^2 ~ " = " ~ .(val)),cex=1.5)
          }
        } 
        par(mai=rep(0,4))
        box(lty=3,col='grey70')
        par(mai=mai)
      }
    }
  }
  options(oop)
  par(opar)    
}

#' \name{pbi_dpairs.legend}
#' \alias{pbi$dpairs.legend}
#' \alias{pbi_dpairs.legend}
#' \title{Adding legend top or bottom to a _dpairs_ or _pairs_ plot.}
#' \description{
#'     The function _dpairs.legend_ allows the user to place a legend outside of a 
#'   pairs or dpairs plot.
#' }
#' \usage{pbi_dpairs.legend(labels,col='grey80',pch=15,cex=1)}
#' \arguments{
#'   \item{labels}{
#'     txt labels to be plotted
#'   }
#'   \item{col}{
#'     colors for the plotting characters
#'   }
#'   \item{pch}{
#'     plotting symbol, default 15
#'   }
#'   \item{cex}{
#'     the character expansion for text and plotting characters, default: 1
#'   }
#' }
#' \examples{
#'   data(iris)
#'   par(omi = c(0.8, 0.4,0.8,0.4)) # reserve some space top and bottom
#'   pbi_dpairs(iris,col.box=2:4,col.xy=rep(c(2:4),each=50))
#'   pbi_dpairs.legend(levels(iris$Species),col=2:4)
#'   mtext('Iris Data',side=3,outer=TRUE,cex=2,line=1)
#' }

pbi_dpairs.legend <- function (labels,col='grey80',pch=15,cex=1) {
  opar=par()
  options(warn=-1)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", labels, xpd = TRUE, horiz = TRUE, inset = c(0,0), 
         bty = "n", pch = pch, col = col, cex = cex)
  par(opar)
}

#' \name{pbi_epsilonSquared}
#' \alias{pbi$epsilonSquared}
#' \alias{pbi_epsilonSquared}
#' \title{Calculate the effect size epsilon-squared for variables of a kruskal.test}
#' \description{
#'     Cohen's rule of thumb for interpretation is: 0.01-0.09 small, 0.09-0.25 medium and above 0.25 large effect. You can convert epsilon-squared to a Spearman _r_ by using the square root of epsilon-square.
#' }
#' \usage{pbi_epsilonSquared(x,y)}
#' \arguments{
#'   \item{x}{
#'     vector with numerical values or a linear model or an aov object.
#'   }
#'   \item{y}{
#'     vector with factor values or numerical vector of a second class
#'   }
#' }
#' \value{return numerical value for epsilon-square for given variables.}
#' \examples{
#'   data(iris)
#'   epsilonSquared=pbi_epsilonSquared
#'   epsilonSquared(iris$Sepal.Length,iris$Species)
#'   # two factor example as well for wilcox.test possible
#'   data(ToothGrowth)
#'   epsilonSquared(ToothGrowth$len,as.factor(ToothGrowth$dose))
#'   # close to r-square of spearman!
#'   cor(ToothGrowth$len,ToothGrowth$dose,method="spearman")^2
#' }
#' \seealso{  See also [pbi](#home), [pbi_etaSquared](#etaSquared) }

pbi_epsilonSquared <- function (x,y) {
  if (class(y) %in% c("numeric","integer")) {
    H=unname(kruskal.test(list(x,y))$statistic)
    n=length(x[which(!is.na(x) & !is.na(y))])
  }  else {
    H=unname(kruskal.test(x ~ y)$statistic)
    n=sum(table(x,y)) # get rid of NA's
    
  }    
  es=H/((n^2-1)/(n+1))
  return(es)
}

#' \name{pbi_etaSquared}
#' \alias{pbi$etaSquared}
#' \alias{pbi_etaSquared}
#' \title{Calculate the effect size eta-squared for an Anova or a linear model. }
#' \description{
#'     Cohen's rule of thumb for interpretation is: 0.01-0.09 small, 0.09-0.25 medium and above 0.25 large effect. You can convert eta-squared to a Pearson r by using the sqrt of eta-square.
#' }
#' \usage{pbi_etaSquared(x,y=NULL)}
#' \arguments{
#'   \item{x}{
#'     vector with numerical values or a linear model or an aov object.
#'   }
#'   \item{y}{
#'     either f factor or NULL if x is given as model.
#'   }
#' }
#' \value{return numerical value for eta-squares for all given variables in the model x or the value for the variable given in y if x is a numerical variable.}
#' \examples{
#'   data(iris)
#'   etaSquared=pbi_etaSquared
#'   etaSquared(iris$Sepal.Length,iris$Species)
#'   etaSquared(lm(iris$Sepal.Length ~ iris$Species))
#'   etaSquared(aov(iris$Sepal.Length ~ iris$Species))
#'   etaSquared(aov(Sepal.Length ~ Species+Sepal.Width,data=iris))
#' }
#' \seealso{  See also [pbi](#home), [pbi_cohensD](#cohensD) }

pbi_etaSquared <- function (x,y=NULL) {
  if (class(x)[1] == "lm") {
    mod=x
    if (length(attr(mod$terms,"dataClasses"))==2) {
      # single factor given
      return(summary(mod)$r.squared)
    } else {
      class(x)="aov"
      return(pbi_etaSquared(x))
    }
  } else if (class(x)[1] == "aov") {
    mod=x
    ss=sum(summary(mod)[[1]][,2])
    sq=summary(mod)[[1]][,2]/ss
    names(sq)=rownames(summary(mod)[[1]])
    sq=sq[1:(length(sq)-1)]
    return(sq)
  } else if (class(x)[1] == "numeric" & class(y)[1] == "factor") {
    mod=aov(x~y)
    return(pbi_etaSquared(mod))
    
  } else {
    stop("Error: wrong call of 'etaSquared'! Call either 'etaSquared(num,factor)' or with 'etaSquared(lm(num~factor))'!")
  }
}

#' \name{pbi_file.head}
#' \alias{pbi$file.head}
#' \alias{pbi_file.head}
#' \title{Displays the first n lines of a file to the terminal}
#' \description{
#'     Displays the first n lines of a file to the terminal
#' }
#' \usage{pbi_file.head(filename,n=6)}
#' \arguments{
#'   \item{filename}{
#'     filename of a text file
#'   }
#'   \item{n}{
#'     number of first lines to display, default: 6
#'   }
#' }
#' \examples{
#'   pbi_file.head(system.file(
#'   "files/pepinfo-spike-sars2.txt",package="pbi"))
#' }
#' \seealso{  See also [sbi$file.cat](#file.cat) }

pbi_file.head = function (filename,n=6) {
  if (!file.exists(filename)) {
    stop(paste('Error! File',filename,'does not exist!'))
  }
  fin=file(filename,'r')
  res=readLines(fin,n=n)
  close(fin)
  return(res)
}

#' \name{pbi_grect}
#' \alias{pbi$grect}
#' \alias{pbi_grect}
#' \title{Adds colored background and grid lines to an existing plot}
#' \description{The function can be used to add a background color to existing plots.
#'   In case the plotting was already done you should add a very transparent color, 
#'   using for instance 33 as the last two digits for a RGB color.
#' }
#' \usage{ pbi_grect(col='#c0c0c033',grid=FALSE) }
#' \arguments{
#'   \item{ col }{
#'     - background color for the plot, default: '#c0c0c033'
#'   }
#'   \item{ grid }{
#'     - should a grid been drawn, default: FALSE
#'   }
#' }
#' \examples{par(mfrow=c(1,2),mai=rep(0.8,4))
#'   data(iris)
#'   pbi_xyplot(iris[,1:2],col=as.numeric(iris$Species)+1)
#'   pbi_grect()
#'   plot(iris[,1:2],type="n")
#'   pbi_grect(col="#ffe0e0",grid=TRUE)
#'   points(iris[,1:2],col=as.numeric(iris$Species)+1,pch=19)
#' }
#' \seealso{  See also: [pbi](#pbi), [pbi_xyplot](#xyplot) }

pbi_grect = function (col="#c0c0c033",grid=FALSE) {
  bods=par("usr")
  rect(bods[1],bods[3],bods[2],bods[4],col=col)
  if (grid) {
    grid (NULL,NULL, lty = 3, col = "grey30")
  }
  
}

#' \name{pbi_impute}
#' \alias{pbi$impute}
#' \alias{pbi_impute}
#' \title{ Missing value imputation using mean, rpart or knn methods.
#' }
#' \description{
#'     Replaces missing values with a reasonable guess by different imputation methods such as 
#'   the simple and not recommended methods mean and median, where NA's are replaced with the 
#'   mean or median for this variable or the more recommended methods using rpart decision trees
#'   or knn using a correlation distance based k-nearest neighbor approach. 
#'   The rpart method can be as well used to replace missing values for categorical variables.
#'   In case of median and mean imputations for categorical variables the modus is used, 
#'   so missing values are replaced with the most often category. This is rarely reasonable.
#' }
#' \usage{pbi_impute(x,method="mean", k=5, cor.method="spearman")}
#' \arguments{
#'   \item{x}{either a matrix or data frame}
#'   \item{method}{ the method used for replacing missing values, either mean, 
#'      median, rpart or knn, default mean (not recommended).
#'   }
#'   \item{k}{in case of knn imputation the number of neighbors, default: 5}
#'   \item{cor.method}{in case of distance determination the correlation method, default: 'spearman'}
#' }
#' \value{return depending on the input either a data frame or matrix with NA's replaced by imputed values.}
#' \examples{
#'   data(iris)
#'   ir=as.matrix(iris[,1:4])
#'   idx=sample(1:length(ir),50)
#'   ir[idx]=NA
#'   summary(ir)
#'   irc=as.matrix(ir)
#'   ir=pbi_impute(ir)
#'   summary(ir)
#'   ir=iris
#'   ir[1,3]=NA; ir[2,4]=NA; ir[c(1,3),5]=NA
#'   head(ir,3)
#'   head(pbi_impute(ir,method="rpart"),3)
#'   irmea=pbi_impute(irc,method="mean")
#'   irknn=pbi_impute(irc,method="knn")
#'   irrpt=pbi_impute(irc,method="rpart")
#'   cor(as.matrix(iris[,1:4])[is.na(irc)],irmea[is.na(irc)])  
#'   cor(as.matrix(iris[,1:4])[is.na(irc)],irknn[is.na(irc)])  
#'   cor(as.matrix(iris[,1:4])[is.na(irc)],irrpt[is.na(irc)])  
#' }
#' \seealso{  See also: [pbi](#pbi) }

pbi_impute <- function (x,method="mean",k=5,cor.method="spearman")  {   
  if (method %in% c("mean","median")) {
    for (i in 1:ncol(x)) {
      # integer is as well numeric :) so 
      if (is.numeric(x[,i])) {
        idx=which(is.na(x[,i]))
        if (method == "mean") {
          x[idx,i]=mean(x[,i],na.rm=TRUE)
        } else if (method == "median") {
          x[idx,i]=median(x[,i],na.rm=TRUE)
        }  
      } else {
        # TODO: modus (?)
        warning(paste("Only numerical columns can be imputed with mean and median! Column",colnames(x)[i], "is however non-numeric!"))
      }
    }
  } else if (method == "rpart") {
    # TODO: refinement for many variables, 
    # take only variables with high absolute correlation
    # into account if more than 10 variables take top 10
    for (i in 1:ncol(x)) {
      idx = which(!is.na(x[,i]))
      if (length(idx) == nrow(x)) {
        next
      }
      if (is.factor(x[,i])) { 
        model=rpart(formula(paste(colnames(x)[i],"~.")), 
                    data=as.data.frame(x[idx,]),
                    method="class")
        x2 = predict(model,newdata=as.data.frame(x[-idx,]),
                     type="class")
      } else {
        model=rpart(formula(paste(colnames(x)[i],"~.")), 
                    data=as.data.frame(x[idx,]))
        x2 = predict(model,newdata=as.data.frame(x[-idx,]))
      }
      
      x[-idx,i]=x2
    }
  } else if (method == "knn") {
    if (ncol(x) < 4) {
      stop("knn needs at least 4 variables / columns")
    }
    data.imp=x
    #D=as.matrix(1-((cor(t(data.imp),use="pairwise.complete.obs")+1)/2))
    D=as.matrix(dist(scale(data.imp)))
    for (i in 1:ncol(x)) {
      idx=which(is.na(x[,i]))
      idxd=which(!is.na(x[,i]))
      for (j in idx) {
        idxo=order(D[j,])
        idxo=intersect(idxo,idxd)
        mn=mean(x[idxo[1:k],i])
        data.imp[j,i]=mn
      }
    }
    return(data.imp)
  } else {
    stop("Unknown method, choose either mean, median, knn or rpart")
  } 
  return(x)
}

#' \name{pbi_lmPlot}
#' \alias{pbi$lmPlot}
#' \alias{pbi_lmPlot}
#' \title{ Plot a xy-plot with linear model fits and the confidence lines. }
#' \description{
#'     The plot visualizes the linear fit description between numerical variables,
#'   together with the confidence limits for the predictions as well as for
#'   the linear model. The code is based on a tutorial by Conrad Halling (2006). 
#'   [https://web.archive.org/web/20180415155316/http://sphaerula.com/legacy/R/linearRegression.html](https://web.archive.org/web/20180415155316/http://sphaerula.com/legacy/R/linearRegression.html)
#' }
#' \usage{ pbi_lmPlot(x,y,col="blue",pch=19,col.lm="red",grid=TRUE,...) }
#' \arguments{
#'   \item{ x }{
#'     vector with numerical variables
#'   }
#'   \item{ y }{
#'     vector with numerical variables
#'   }
#'   \item{ col }{
#'     color for points and confidence lines for predictions, default: blue
#'   }
#'   \item{ pch }{
#'     plotting character, default: 19
#'   }
#'   \item{col.lm}{
#'     color for the linear model and the confidence lines, default: red
#'   }
#'   \item{ grid }{
#'     should a plot grid be drawn, default: TRUE
#'   }
#'   \item{\ldots}{other arguments delegated to the standard plot function}
#' }
#' \examples{
#' c20.22=c(
#' 17.9, 18.3, 18.3, 18.4, 18.4, 20.2, 20.3, 21.8, 21.9, 
#' 22.1, 23.1, 24.2, 24.4)
#' ins.sens=c(
#' 250, 220, 145, 115, 230, 200, 330, 400, 370, 260, 270, 
#'  530, 375)
#' 
#' pbi_lmPlot(x=c20.22,y=ins.sens,
#'    xlab='\%C20-22 Fatty Acids',ylim=c(0,600),
#'    xlim=c(17,25),main='best fit',
#'    ylab='Insuline Sensitivity Index (mg/m^2/min)')
#' 
#' legend('bottomright',c('best fit','fit 95\% CI',
#' 'prediction 95\% CI'),lty=c(1,2,2),
#' col=c('red','red','blue'))
#' }
#' \seealso{  See also: [pbi](#pbi), [pbi_modelQuality](#modelQuality) }

pbi_lmPlot = function (x,y, col="blue",pch=19,col.lm="red",grid=TRUE,...) {
  df <- data.frame(x=x,y=y)
  plot(y ~ x, data=df, pch=pch, col=col,...)
  if (grid) {
    grid (NULL,NULL, lty = 3, col = "grey30")
  }
  mod <- lm( y ~ x, data=df)
  new.x.df = data.frame(
    x = seq(from   = range(df$x)[1],
            to     = range(df$x)[2],
            length = 11 ))
  lty=list(upr=2,lwr=2,fit=1)
  for (lim in c("upr","lwr","fit")) {
    lines(x   = new.x.df$x,
          y   = predict(mod, new.x.df, 
                        interval = "confidence" )[,lim],
          col = col.lm,lty=lty[[lim]],lwd=2 )
  }
  for (lim in c("upr","lwr")) {
    lines(
      x   = new.x.df$x,
      y   = predict( mod, new.x.df, 
                     interval = "prediction" )[ , lim ],
      col = col ,lty=2,lwd=2)
  }
}

#' \name{pbi_mi}
#' \alias{pbi$mi}
#' \alias{pbi_mi}
#' \title{ Return the mutual information for two vectors or a binned table. }
#' \description{
#'     Returns: mutual information value, or matrix of all pairwise values in case input  is matrix or data frame
#' }
#' \usage{ pbi_mi(x,y=NULL, breaks=4) }
#' \arguments{
#'   \item{ x }{
#'     either a binned table, a numerical vector, a data frame or matrix
#'   }
#'   \item{ y }{
#'     a numerical vector if x is not a binned table, matrix or data frame
#'   }
#'   \item{ breaks }{
#'     number of breaks to create a binned table if x and y are numerical vectors, default: 4
#'   }
#' }
#' \value{return mutual information value or matrix of all pairwise values in case input is matrix or data frame.}
#' \examples{
#' rn1=rnorm(100,mean=10,sd=1);
#' rn2=rn1+0.5*rnorm(100)
#' cor(rn1,rn2) # high
#' cor(rn1,sample(rn2)) #random 
#' pbi_mi(rn1,rn2) # high 
#' pbi_mi(rn1,sample(rn2)) #random
#' tab=table(cut(rn1,breaks=4),cut(rn2,breaks=4))
#' pbi_mi(tab)
#' pbi_mi(rn1,rn2,breaks=4)
#' pbi_mi(rn1,rn2,breaks=7)
#' data(iris)
#' round(pbi_mi(iris[,1:4]),2)
#' mii=pbi_mi(iris[,1:4])
#' round(mii/diag(mii),2)
#' round(cor(iris[,1:4]),2)
#' }

pbi_mi = function (x,y=NULL,breaks=4) {
  if (is.data.frame(x) | is.matrix(x)) {
    M=matrix(0,nrow=ncol(x),ncol=ncol(x))
    colnames(M)=rownames(M)=colnames(x)
    for (i in 1:ncol(x)) {
      for (j in 1:ncol(x)) {
        if (i>=j) {
          M[i,j]=M[j,i]=pbi_mi(x[,i],x[,j],breaks=breaks)
        }
      }
    }
    return(M)
  } else {
    if (!(class(x)=='table')) {
      x=table(cut(x,breaks=breaks),cut(y,breaks=breaks))        
    }
    f1=x/sum(x)
    fx=rowSums(f1)
    fy=colSums(f1)
    fn=fx %o% fy
    f2=fn/sum(fn)
    LR = ifelse(f1>0,log(f1/f2),0)
    MI = sum(f1*LR)
    return(MI)
  }
}

#' \name{pbi_modelQuality}
#' \alias{pbi$modelQuality}
#' \alias{pbi_modelQuality}
#' \title{ Quality measures for determining model qualitiy. }
#' \description{
#'   The functions returns a few basic model quality measures 
#'   such as RMSE, MSE and MAE in case the data in x and y are numeric, 
#'   or Accuracy, Sensitivity, Specificity and Balanced Classification rate (BCR)
#'   in case the data are factors.
#' }
#' \usage{ pbi_modelQuality(x,y) }
#' \arguments{
#'   \item{ x }{
#'     numercial values or factor values representing the true values.
#'   }
#'   \item{ y }{
#'     numerical values or factor values representing the predicted values.
#'   }
#' }
#' \value{return list with the components MAE, MSE, RMSE for numerical values or ACC, SENS, SPEC and BCR for factors (classes).}
#' \examples{
#' library(MASS)
#' data(swiss)
#' # regression measures
#' mod=lm(bwt~lwt,data=birthwt)
#' summary(mod)$r.squared
#' pred=predict(mod,newdata=birthwt)
#' unlist(pbi_modelQuality(birthwt$bwt,pred))
#' mod=lm(bwt~lwt+age,data=birthwt)
#' summary(mod)$r.squared
#' unlist(pbi_modelQuality(birthwt$bwt,predict(mod,newdata=birthwt)))
#' 
#' # classification measures
#' library(rpart)
#' # simple two class problem
#' sex=as.factor(c(rep("M",50),rep("F",50)))
#' height=c(rnorm(50,mean=180,sd=7),rnorm(50,mean=178,sd=6))
#' mod=rpart(sex~height)
#' table(sex,pred=predict(mod,newdata=data.frame(height=height),type="class"))
#' unlist(pbi_modelQuality(sex,
#'        predict(mod,newdata=data.frame(height=height),type="class")))
#' # simple two class problem but unbalanced
#' sex=as.factor(c(rep("M",30),rep("F",70))) # rare males
#' height=c(rnorm(30,mean=180,sd=7),rnorm(70,mean=178,sd=6))
#' mod=rpart(sex~height)
#' table(sex,pred=predict(mod,newdata=data.frame(height=height),type="class"))
#' unlist(pbi_modelQuality(sex,
#'        predict(mod,newdata=data.frame(height=height),type="class")))
#' # difficult three class problem
#' table(iris$Species)
#' mod=rpart(Species~.,data=iris)
#' pred=predict(mod,newdata=iris,type="class")
#' table(pred,iris$Species)
#' # convert to two classes (setosa)
#' table(c('s','v','v')[as.numeric(iris$Species)])
#' unlist(pbi_modelQuality(
#'    as.factor(c('s','v','v')[as.numeric(iris$Species)]),
#'    as.factor(c("s","v","v")[as.numeric(pred)])))
#' # convert to two classes (versicolor)
#' unlist(pbi_modelQuality(
#'    as.factor(c('x','v','x')[as.numeric(iris$Species)]),
#'    as.factor(c("x","v","x")[as.numeric(pred)])))
#' # convert to two classes (virginica)
#' unlist(pbi_modelQuality(
#'    as.factor(c('x','x','v')[as.numeric(iris$Species)]),
#'    as.factor(c("x","x","v")[as.numeric(pred)])))
#' }
#' \seealso{  See also: [pbi](#pbi), [pbi_lmPlot](#lmPlot) }

pbi_modelQuality = function (x,y) {
  mse = function (x,p) { # outlier sensitive
    sum=sum((x-p)^2,na.rm=TRUE)
    N=length(which(!is.na(x+p)))
    return(sum/N)
  }
  mae = function (x,p) { #easier to interpete
    sum=sum(abs(x-p),na.rm=TRUE)
    N=length(which(!is.na(x+p)))
    return(sum/N)
  }   
  if (!length(x) == length(y)) {
    stop("Error: Both vectors x and y must have the same length!")
  }
  if (is.numeric(x) & !is.numeric(y)) {
    stop("Error: Both vectors x and y must have the same type!")
  } else if (is.numeric(x)) {
    return(list(MAE=mae(x,y),MSE=mse(x,y),RMSE=sqrt(mse(x,y)),r.squared=cor(x,y)^2))
  } else if (is.factor(x) & !is.factor(y)) {
    y=factor(y,levels=levels(y))
    return(pbi_modelQuality(x,y))
  } else if (is.factor(x)) {
    tab=table(y,x)   
    TP=tab[1,1]; FP=tab[1,2]; FN=tab[2,1]; TN=tab[2,2]
    ACC=(TP+TN)/sum(tab)
    SENS=TP/(TP+FN)
    SPEC=TN/(TN+FP)
    BCR=(SENS+SPEC)/2
    return(list(ACC=ACC,SENS=SENS,SPEC=SPEC,BCR=BCR))
  }
}

#' \name{pbi_modus}
#' \alias{pbi$modus}
#' \alias{pbi_modus}
#' \title{ Return the most often level in a categorical variable. }
#' \description{ Simple implementation of returning the most most often appearinglevel(s) of a categorical variable. }
#' \usage{ pbi_modus(x) }
#' \arguments{
#'   \item{ x }{
#'     a vector with elements of class factor
#'   }
#' }
#' \value{return most often apparent level in the categorical variable}
#' \examples{
#'   pbi_modus(c('A','A','B','C'))
#'   pbi_modus(c('A','A','B','B','C'))
#' }

pbi_modus = function (x) {
    cat=x
  tab=table(cat)
  idx=which(max(tab)==tab)
  return(names(tab)[idx])
}

#' \name{pbi_msa2pwm}
#' \alias{pbi$msa2pwm}
#' \alias{pbi_msa2pwm}
#' \title{Calculates PFM, PPM and PWM for the given alignment file.}
#' \description{
#'     Implementation is partially derived from Dave Tang, see: https://davetang.org/muse/2013/10/01/position-weight-matrix/
#' }
#' \usage{pbi_msa2pwm(filename)}
#' \arguments{
#'   \item{filename}{
#'     filename of a file wither with two columns (second column with aligned sequences, or a single column file with the first column having the sequences.
#'   }
#' }
#' \value{return list with the following components:
#'   - PFM - Position frequency matrix
#'   - PPM - Position probability matrix
#'   - PWM - Position weight matrix
#'   - PWMPC - Position weight matrix with added pseudocounts (untested)}
#' \examples{
#' # data example: https://en.wikipedia.org/wiki/Position_weight_matrix
#'   cat("GAGGTAAAC
#' TCCGTAAGT
#' CAGGTTGGA
#' ACAGTCAGT
#' TAGGTCATT
#' TAGGTACTG
#' ATGGTAACT
#' CAGGTATAC
#' TGTGTGAGT
#' AAGGTAAGT
#' ",file="test.clu")
#' lapply(pbi_msa2pwm("test.clu"),round,2)
#' }

pbi_msa2pwm <- function (filename) {
  ## https://davetang.org/muse/2013/10/01/position-weight-matrix/
  pwm <- function(freq,total,bg=0.25){
    #using the formulae above
    p <- (freq + (sqrt(total) * bg)) / (total + (1/bg * (sqrt(total) * bg)))
    return(log2(p/bg))
  }
  fin  = file(filename, "r")
  l=c()
  ncol=0
  while(length((line = readLines(fin,n=1)))>0) {
    if (grepl("^[A-Za-z0-9][^\\s]+ +[A-Za-z]+$",line)) {
      seq=gsub("^[^ ]+ +([A-Za-z]+)$","\\1",line)
      l = c(l,strsplit(toupper(seq),"")[[1]])
      ncol=nchar(seq)
    } else  if (grepl("^[A-Za-z]+ *$",line)) {
      seq=gsub("^([^ ]+) *$","\\1",line)
      l = c(l,strsplit(toupper(seq),"")[[1]])
      ncol=nchar(seq)
    }
  }
  close(fin)
  M=matrix(l,ncol=ncol,byrow=TRUE)
  PPM=matrix(0,nrow=26,ncol=ncol)
  PFM=matrix(0,nrow=26,ncol=ncol)
  rownames(PPM)=rownames(PFM)=LETTERS
  for (i in 1:ncol) {
    T=table(M[,i])
    PFM[names(T),i]=as.numeric(T)
    T=prop.table(T)
    PPM[names(T),i]=as.numeric(T)
  }
  idx=which(apply(PPM,1,sum)>0)
  bg=0.25
  if (length(idx)>4) {
    bg=0.05
  }
  PWM=pwm(PFM[idx,],sum(PFM[idx,1]),bg=bg)
  return(list(PFM=PFM[idx,],PPM=PPM[idx,],PWM=log2(PPM[idx,]/bg),PWMPC=PWM))
}

#' \name{pbi_mw}
#' \alias{pbi$mw}
#' \alias{pbi_mw}
#' \title{ Determine the molecular weight for a given sequence (sequences). }
#' \description{Calculates and returns the molecular weight for a given amino acid sequence.}
#' \usage{ pbi_mw(seq)}
#' \arguments{
#'   \item{ seq }{
#'     an amino acid sequence string
#'   }
#' }
#' \value{return the protein molecular weight.}
#' \examples{
#'   pbi_mw("AACTLIS")
#'   pbi_mw("A")
#' }
pbi_mw <- function (seq) {
  tab=read.table(text='
  Alanine	Ala	A	89
  Arginine	Arg	R	174
  Asparagine	Asn	N	132
  Aspartic-acid	Asp	D	133
  Cysteine	Cys	C	121
  Glutamine	Gln	Q	146
  Glutamic-acid	Glu	E	147
  Glycine	Gly	G	75
  Histidine	His	H	155
  Isoleucine	Ile	I	131
  Leucine	Leu	L	131
  Lysine	Lys	K	146
  Methionine	Met	M	149
  Phenylalanine	Phe	F	165
  Proline	Pro	P	115
  Serine	Ser	S	105
  Threonine	Thr	T	119
  Tryptophan	Trp	W	204
  Tyrosine	Tyr	Y	181
  Valine	Val	V	117',row.names=3)
  colnames(tab)=c('Name','TLC','Da')
  mw=0
  for (aa in strsplit(seq,'')[[1]]) {
    mw=mw+tab[aa,'Da']
  }   
  return(mw)
}

#' \name{pbi_package.deps}
#' \alias{pbi$package.deps}
#' \alias{pbi_package.deps}
#' \title{ Return the packages which are required by the given package name. }
#' \description{The function helps in identifying the nested package dependencies for a given package.}
#' \usage{ pbi_package.deps(pkgName,mode='all') }
#' \arguments{
#'   \item{ pkgName }{
#'     an package name given as text string.
#'   }
#'   \item{ mode }{
#'     which package names to return, the following modes are available:
#'     'all' - all required packages,  
#'     'install' - not yet installed packages, 
#'     'nonbase' - packages not in the standard R installation
#'   }
#' }
#' \value{return list of required packages.}
#' \examples{
#'   \dontrun{
#'     pbi_package.deps('igraph',mode='nonbase')
#'     pbi_package.deps('igraph',mode='all')
#'  }
#' }

pbi_package.deps <- function(pkgName,mode='all')  {
  x=pkgName
  if (!interactive()) {
    r <- getOption("repos");
    r["CRAN"] <- "https://www.freestatistics.org/cran/"
    #"https://lib.ugent.be/CRAN/" 
    options(repos=r) 
  }
  deps=tools::package_dependencies(x,recursive=TRUE)[[1]]
  if (mode == 'install') {
    idx = which(
      !(deps %in% rownames(installed.packages())))
    return(deps[idx])
  } else if (mode == 'nonbase') {
    ipacks=installed.packages()
    bpacks=ipacks[ipacks[,'Priority'] %in% 
                    c('base','recommended'),]
    rnms=setdiff(rownames(ipacks),rownames(bpacks))
    return(intersect(deps,rnms))
  } else if (mode == 'all') {
    return(deps)
  } else {
    stop('mode must be either `all`, `install` or `nonbase`!')
  }
  
}

#' \name{pbi_pastel}
#' \alias{pbi$pastel}
#' \alias{pbi_pastel}
#' \title{ Create up to 20 pastel colors }
#' \description{
#'   This is an alternative color creation function for R versions before 3.6 where 
#'   the function `hcl.colors` is not available.
#' }
#' \usage{ pbi_pastel(n) }
#' \arguments{
#'   \item{ n }{
#'     number of colors requested, must be within 2 and 20
#'   }
#' }
#' \value{return vector of colors in RGB codes of requested length 'n'.}
#' \examples{
#' pbi_pastel(4)
#' par(mai=c(0.2,0.2,0.2,0.1))
#' plot(1:20,col=pbi_pastel(20),cex=3,pch=15)
#' }

pbi_pastel <- function (n) {
  if(n > 20 |  n < 1) {
    stop("only between 1 and 20 colors can be given" ) 
  }
  pcols= c("#FFC5D0","#FDC8C3","#F6CBB7","#EDD0AE","#E2D4A8","#D4D8A7","#C5DCAB","#B6DFB4","#A8E1BF",
           "#9EE2CB", "#99E2D8","#9BE0E5","#A4DDEF","#B3D9F7","#C4D5FB","#D5D0FC","#E4CBF9","#F0C7F2",
           "#F9C5E9", "#FEC4DD")
  idx=seq(1,20,by=floor(20/n))
  return(pcols[idx])
}
#' \name{pbi_pca.biplot}
#' \alias{pbi$pca.biplot}
#' \alias{pbi_pca.biplot}
#' \title{Improved biplot for pca objects.}
#' \description{
#'   The function _pca.biplot_ provides an improved biblot for
#'   visualizing the pairwise scores of individual principal components of 
#'   an object created using the function _prcomp_. In contrast to the default 
#'   biplot function  this plot visualizes the data as points and not row numbers,
#'   it allows to display groups using color codes and distribution ellipses.
#' }
#' \usage{pbi_pca.biplot(pca,pcs=c("PC1","PC2"),pch=19,col='black',text=NULL,arrows=TRUE,arrow.fac=1,arrow.n=-1,ellipse=FALSE,ell.fill=FALSE,xlab=NULL,ylab=NULL,grid=TRUE,scale=NULL,...)}
#' \arguments{
#'   \item{pca}{
#'     pca object of class _prcomp_, created using the function _prcomp_.
#'   }
#'   \item{pcs}{
#'     the components to plot, default: c('PC1','PC2')
#'   }
#'   \item{pch}{
#'     plotting character, default: 19
#'   }
#'   \item{col}{
#'     plotting color, default: black
#'   }
#'   \item{text}{
#'     optional text labels, if not given  pch is used, default: NULL
#'   }
#'   \item{arrows}{
#'     should loading arrows be displayed, default: TRUE
#'   }
#'   \item{arrow.fac}{
#'     scaling factor for arrow length, default: 1
#'   }
#'   \item{arrow.n}{
#'     how many arows per PC to display, default: -1 (all)
#'   }
#'   \item{ellipse}{
#'     should 85 and 95 confidence intervals for the chisq distribution be shown. If this is shown colors for each group using the col argument must be given, default: FALSE
#'   }
#'   \item{ell.fill}{
#'     should a filled 85 percent confidence interval be shown, colors will be used from the plotting color with opacity, default: FALSE
#'   }
#'   \item{xlab}{
#'     custom xlab, if not given the PC name with variance in percent is shown, default: NULL
#'   }
#'   \item{ylab}{
#'     custom ylab, if not given the PC name with variance in percent is shown, default: NULL
#'   }
#'   \item{grid}{
#'     should a plot grid be drawn, default: TRUE
#'   }
#'   \item{scale}{
#'     function to scale to coordinate values sich as asinh, default: NULL
#'   }
#'   \item{...}{
#'     additional arguments delegated to the standard plot function
#'   }
#' }
#' \examples{
#'   par(mfrow=c(1,2),mai=c(0.8,0.8,0.6,0.2))
#'   data(iris)
#'   pci=prcomp(iris[,1:4],scale=TRUE)
#'   pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=TRUE,ell.fill=TRUE,
#'       arrow.fac=2.3,arrows=TRUE)
#'   legend('topright',pch=19,col=2:4,levels(iris$Species))
#'   # just a score plot
#'   pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=TRUE,ell.fill=TRUE,
#'       arrow.fac=2.3,arrows=FALSE)
#'   pbi_pca.biplot(pci,col=rep(2:4,each=50),ellipse=FALSE,arrow.fac=2.3,arrows=FALSE)
#'   pbi_pca.biplot(pci,pcs=c('PC1','PC3'),col=rep(2:4,each=50),
#'          ellipse=FALSE,arrow.fac=2.3,arrows=FALSE)
#'   data(swiss)
#'   col=c(2,4)[as.numeric(cut(swiss$Catholic,breaks=c(0,20,100)))]
#'   pcs=prcomp(swiss,scale=TRUE)
#'   pbi_pca.biplot(pcs,arrow.fac=2,grid=FALSE,
#'         col=col,ellipse=TRUE,text=substr(rownames(swiss),1,3))
#'   pbi_pca.biplot(pcs,arrow.fac=2,grid=TRUE, 
#'        col=col,text=substr(rownames(swiss),1,3),scale=asinh,arrows=FALSE)
#' }

pbi_pca.biplot = function (pca,pcs=c("PC1","PC2"),
                             pch=19,col='black',text=NULL,
                             arrows=TRUE,arrow.fac=1,arrow.n=-1,
                             ellipse=FALSE,ell.fill=FALSE,xlab=NULL,ylab=NULL,
                             grid=TRUE,scale=NULL,...) {
  if (missing("xlab")) {
    xlab=paste(pcs[1]," (", round(summary(pca)$importance[2,pcs[1]]*100,1),"%)",sep="")
  }
  if (missing("ylab")) {
    ylab=paste(pcs[2]," (", round(summary(pca)$importance[2,pcs[2]]*100,1),"%)",sep="")
  }
  if (class(scale)=="function") {
    x=scale(pca$x[,pcs[1]])
    y=scale(pca$x[,pcs[2]])
    xlab=paste(substitute(scale),xlab)
    ylab=paste(substitute(scale),ylab)        
    #plot(x,y,pch=pch,col=col,type="n",xlab=xlab,ylab=ylab,...)
  } else {
    x=pca$x[,pcs[1]]
    y=pca$x[,pcs[2]]
  }
  pcdf=data.frame(x=x,y=y)
  plot(x,y,pch=pch,col=col,type="n",xlab=xlab,ylab=ylab,...)
  colnames(pcdf)=c(pcs[1],pcs[2])
  abline(h=0,lty=2)
  abline(v=0,lty=2)    
  if (ellipse) {
    if (length(col)!= nrow(pca$x)) {
      stop("colors must have same length as data points")
    }
    ell.col=col
    i=1
    for (cl in names(table(ell.col))) {
      C=cov(pcdf[ell.col==cl,c(pcs[1],pcs[2])])    # Covarianz-Matrix C bestimmen
      d85=qchisq(0.85, df = 2)     # 85% - Faktor , um die Ellipse zu skalieren
      M=colMeans(pcdf[ell.col==cl,c(pcs[1],pcs[2])]) #   Mittelwerte (Zentrum) des Clusters
      el=cluster::ellipsoidPoints(C, d85, loc=M)  # Ellipsen-Punkte aus C und M berechnen
      if (ell.fill) {
        colfill=paste(rgb(t(col2rgb(cl))/255),"33",sep="")
        polygon(el,col=colfill,border=NA)
        i=i+1
        next
      }
      lines(el,col=cl,lwd=1.5,lty=2)    #  Ellipse als geschlossene Linies zeichnen
      d95=qchisq(0.95, df = 2)     # 85% - Faktor , um die Ellipse zu skalieren
      el=cluster::ellipsoidPoints(C, d95, loc=M)  # Ellipsen-Punkte aus C und M berechnen
      lines(el,col=cl,lwd=1.5,lty=1)    #  Ellipse als geschlossene Linies zeichnen                        
      
    }
  }
  if (class(text[1])!="NULL") {
    if (length(text)==1) {
      text=rownames(pca$x)
      
    } 
    text(pcdf[,pcs[1]],pcdf[,pcs[2]],text,col=col,...)
  } else {
    points(pcdf[,pcs[1]],pcdf[,pcs[2]],pch=pch,col=col,...)
  }
  if (arrows) {
    loadings=pca$rotation
    if (arrow.n == -1) {
      idx=1:nrow(loadings)
    } else {
      load1=loadings[order(abs(loadings[,pcs[[1]]])),]
      load2=loadings[order(abs(loadings[,pcs[[2]]])),]
      idx=which(rownames(loadings) %in% c(rownames(load1)[1:arrow.n],
                                          rownames(load2)[1:arrow.n]))
    }
    arrows(0,0,loadings[idx,pcs[1]]*arrow.fac,loadings[idx,pcs[2]]*arrow.fac,
           length=0.1,angle=20,col='black')
    text(loadings[idx,pcs[1]]*arrow.fac*1.2,loadings[idx,pcs[2]]*arrow.fac*1.2,
         rownames(loadings)[idx],col='black',font=2)
  }
  if (grid) {
    grid (NULL,NULL, lty = 3, col = "grey30")
  }
  
}

#' \name{pbi_pca.corplot}
#' \alias{pbi_pca.corplot}
#' \alias{pbi$pca.corplot}
#' \title{PCA correlation plot to show association between PCs and variables.}
#' \description{
#'   The function provides a PCA correlation plot to show associations 
#'   between PCs and variables. The closer a variable to the PC coordinate 
#'   the higher the correlation, the more away from the center of the 
#'   coordinate system, the higher the impact of the variable on this PC.
#'   You can think about the corplot as a biplot ommiting the samples and 
#'   the arrows.
#' }
#' \usage{pbi_pca.corplot(pca,pcs=c("PC1","PC2"), main="Correlation plot",cex=NULL,nvar=64,...)}
#' \arguments{
#'   \item{pca}{
#'     pca object which was created using the function _prcomp_.
#'   }
#'   \item{pcs}{
#'     vector of two PCs to be plotted against each other, default: c('PC1','PC2')
#'   }
#'   \item{main}{
#'     title of the plot, default: 'Correlation plot'
#'   }
#'   \item{cex}{
#'     character expansion for the samples, default: NULL (automatic calculation)
#'   }
#'   \item{nvar}{
#'     number of variables which will be displayed, for both components the variables which contributes mostly to the variances will be used, default: 64.
#'   }
#'   \item{...}{
#'     remaining arguments are delegated to the standard plot function
#'   }
#' }
#' \examples{
#' par(mfrow=c(1,2),mai=c(0.7,0.7,0.5,0.1))
#' library(cluster)
#' data(votes.repub)
#' pca=prcomp(t(na.omit(votes.repub)))
#' pbi_pca.corplot(pca)
#' data(swiss)
#' pca=prcomp(swiss,scale.=TRUE)
#' pbi_pca.corplot(pca)
#' }

pbi_pca.corplot = function (pca,pcs=c("PC1","PC2"), main="Correlation plot",cex=NULL,nvar=64,...) {
  if (is.logical(pca$scale)) {
    df=t(t(pca$x %*% t(pca$rotation)) + pca$center)
  } else {
    df=t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)  
  }
  xlab=paste(pcs[1]," (", round(summary(pca)$importance[2,pcs[1]]*100,1),"%)",sep="")
  ylab=paste(pcs[2]," (", round(summary(pca)$importance[2,pcs[2]]*100,1),"%)",sep="")    
  getMainLoadings = function (pca,PC="PC1",n=10) {
    return(rownames(head(pca$rotation[rev(order(abs(pca$rotation[,PC]))),],n)))
  }
  #this=dpca
  plot(1,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),pch="",
       xlab=xlab,
       ylab=ylab,
       asp=1,main=main,...)
  # library plotrix       
  draw.circle =  function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
                           lwd = 1) 
  {
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
      col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
      xv <- cos(angles) * radius[circle] + x
      yv <- sin(angles) * radius[circle] * ymult + y
      polygon(xv, yv, border = border, col = col[circle], lty = lty, 
              lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
  }
  draw.circle(0,0,1,lty=2)
  lines(c(0,0),c(1.1,-1.1),lwd=2)
  lines(c(1.1,-1.1),c(0,0),lwd=2)
  x=1
  if (length(colnames(df)>nvar)) {
    cnames=getMainLoadings(pca,pcs[1],nvar)
    cnames2=getMainLoadings(pca,pcs[2],nvar)
    cnames=unique(cnames,cnames2)
  } else {
    cnames =colnames(df)
  }
  if (is.null(cex)) {
    text.cex=0.4+4/length(cnames)
    if (text.cex > 1.2) {
      text.cex=1.2
    }
  } else {
    text.cex=cex
  }
  mcor1=c()
  mcor2=c()
  for (cname in cnames) {
    mcor1=c(mcor1,cor.test(df[,cname],pca$x[,pcs[1]])$estimate)
    mcor2=c(mcor2,cor.test(df[,cname],pca$x[,pcs[2]])$estimate)
  }
  text(mcor1,mcor2,labels=cnames,cex=text.cex)
}

#' \name{pbi_pca.pairs}
#' \alias{pbi$pca.pairs}
#' \alias{pbi_pca.pairs}
#' \title{Improved pairs plot for pca objects.}
#' \description{
#'   The function _pca.pairs_ provides an improved pairs plot for
#'   visualizing the pairwise scores of the individual components of an analyses 
#'   using the function _prcomp_. In contrast to the default pairs function 
#'   this plot visualizes in the diagonal as well the variances and 
#'   a density line for the component scores.
#' }
#' \usage{pbi_pca.pairs(pca,n=10,groups=NULL,col='black',pch=19,legend=FALSE,...)}
#' \arguments{
#'   \item{pca}{
#'     pca object which was created using the function _prcomp_.
#'   }
#'   \item{n}{
#'     maximal number of components to visualize, default: 10
#'   }
#'   \item{groups}{
#'     vector with classes having the same length than the input matrix for prcomp has rows, default: NULL
#'   }
#'   \item{col}{
#'     colors for the plotting, character, default: 'black'
#'   }
#'   \item{pch}{
#'     plotting, symbol, default: 19
#'   }
#'   \item{legend}{
#'     should the legend be displayed on top, default: FALSE
#'   }
#'   \item{...}{
#'     additional arguments delegated to the standard _pairs_ function
#'   }
#' }
#' \examples{
#'   data(iris)
#'   pci=prcomp(iris[,1:4],scale=TRUE)
#'   pbi_pca.pairs(pci,pch=15,groups=iris[,5],
#'      legend=TRUE,oma=c(5,4,4,4),col=as.numeric(iris[,5])+1)
#' }

N = new.env()
N$.n = 0
pbi_pca.pairs = function (pca,n=10,groups=NULL,
                            col='black',pch=19,legend=FALSE,...) {
  N$.n <- 1
  if (n>ncol(pca$x)) {
    n=ncol(pca$x)
  }
  pst=FALSE
  if (class(groups) != "NULL" & length(col) != length(groups)) {
    coln=length(levels(as.factor(groups)))
    cols=pbi_pastel(coln)
    col=cols[as.numeric(as.factor(as.character(groups)))]
    pst=TRUE
  }
  panel.text = function (x,...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    text(0.5,0.5,
         paste(sprintf("%.1f",summary(pca)$importance[2,  N$.n]*100),"%",
               sep=""),cex=1.5)
    ds=density(pca$x[,N$.n],na.rm=TRUE)
    ds$x=ds$x-min(ds$x)
    ds$x=ds$x/max(ds$x)
    ds$y=(ds$y/max(ds$y)*0.3)
    polygon(ds,col='grey80')
    N$.n = N$.n + 1
  }
  pairs(pca$x[,1:n],diag.panel=panel.text,col=col,pch=pch,...)
  if (legend && class(groups) != "NULL") {
    opar=par()
    options(warn=-1)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    if (pst) {
      leg.cols=pbi_pastel(coln)[as.numeric(as.factor((levels(as.factor(groups)))))]
    } else {
      cols=col
      names(cols)=as.character(groups)
      lcol=cols[unique(names(cols))]
      leg.cols=as.numeric(lcol)
    }
    legend('bottom', levels(as.factor(groups)), xpd = TRUE, 
           horiz = TRUE, inset = c(0,0), 
           bty = "n", pch = pch, col = leg.cols, cex = 1.2)
    
    par(opar)
  }
  
}

#' \name{pbi_pca.plot}
#' \alias{pbi$pca.plot}
#' \alias{pbi_pca.plot}
#' \title{Improved bar or screeplot for pca objects.}
#' \description{
#'   The function _pca.plot_ provides an improved bar- or screeplot for
#'   visualizing the variances of the individual components of an analyses 
#'   using the function _prcomp_. In contrast to the default plot function 
#'   this plot visualize cumulative and individual variances in percent.
#' }
#' \usage{pbi_pca.plot(pca,n=10,type="bar", cex=1.5, legend=TRUE,xlab="Components",ylab="Variance (\%)",pc.col=c("light blue","grey"),...)}
#' \arguments{
#'   \item{pca}{
#'     pca object which was created using the function _prcomp_.
#'   }
#'   \item{n}{
#'     maximal number of components to visualize, default: 10
#'   }
#'   \item{type}{
#'     plotting type either "bar" or "scree", default: "bar"
#'   }
#'   \item{cex}{
#'     character expansion for the legend and the screeplot plotting characters, default: 1.5
#'   }
#'   \item{legend}{
#'     should the legend be displayed on top, default: TRUE
#'   }
#'   \item{xlab}{
#'     
#'   }
#'   \item{ylab}{
#'     
#'   }
#'   \item{pc.col}{
#'     colors for the PC variances, first individual, second color for the cumulative variance, default: c("light blue","grey")
#'   }
#'   \item{...}{
#'     additional arguments delegated to the standard plot function
#'   }
#' }
#' \examples{
#'   data(iris)
#'   par(mfrow=c(1,2))
#'   pcai=prcomp(iris[,1:4],scale=TRUE)
#'   pbi_pca.plot(pcai)
#'   pbi_pca.plot(pcai,type="scree",legend=FALSE)
#' }

pbi_pca.plot = function (pca,n=10,type="bar", cex=1.5, 
                           legend=TRUE,xlab="Components",ylab="Variance (%)",
                           pc.col=c("light blue","grey"),...) {
  if (n>ncol(pca$x)) {
    n=ncol(pca$x)
  }
  if (legend) {
    ylim=c(0,120)
  } else {
    ylim=c(0,105)
  }
  if (type=="bar") {
    barplot(summary(pca)$importance[3,1:n]*100,
            ylim=ylim,col='white',
            xlab=xlab,ylab=ylab,axes=FALSE,...)
  } else {
    plot(summary(pca)$importance[3,1:n]*100,type="b",
         ylim=ylim,cex.axis=1.2,lwd=2,cex=cex,
         xlab=xlab,ylab=ylab,axes=FALSE,
         pch=15,col=pc.col[2],...)
    points(summary(pca)$importance[2,1:n]*100,type="b",cex=cex,
           lwd=2,xlab="", pch=15,col=pc.col[1],...)
    
    axis(1,at=1:n,labels=paste("PC",1:n,sep=""))
  }
  axis(2,at=c(20,40,60,80,100),labels=c(20,40,60,80,100))
  if (type == "bar") {
    barplot(summary(pca)$importance[3,1:n]*100,add=TRUE,col=pc.col[2],axes=FALSE)
    barplot(summary(pca)$importance[2,1:n]*100,add=TRUE,col=pc.col[1],axes=FALSE)        
  }
  abline(h=5,lty=2,lwd=0.5)
  abline(h=10,lty=2,lwd=0.5)
  abline(h=90,lty=2,lwd=0.5)
  abline(h=95,lty=2,lwd=0.5)    
  abline(h=100,lty=1,lwd=0.5)    
  if (legend) {
    legend("topleft",c("Component","Cumulative"),col=pc.col,pch=15,cex=1.5,box.lwd=0,ncol=2)
  }
  box()
}

#' \name{pbi_pca.variances}
#' \alias{pbi$pca.variances}
#' \alias{pbi_pca.variances}
#' \title{Return the absolute variance contributions for each variable to each component.}
#' \description{
#'   The function eturn the absolute variance contributions for each variable to each component.
#'   Every squared loading value for each component and variable  is multiplied
#'   with the component importance. The sum of the returned matrix is therefor 1.
#' }
#' \usage{pbi_pca.variances(pca)}
#' \arguments{
#'   \item{pca}{
#'     a PCA object created with `prcomp`.
#'   }
#' }
#' \value{return matrix with absolute variances for each component and variable.}
#' \examples{
#' pca=prcomp(USArrests,scale=TRUE)
#' round(pbi_pca.variances(pca),2)
#' sum(pbi_pca.variances(pca))
#' }
#' \seealso{  See also [pbi](#pbi), [pbi$pca.varplot](#pca.varplot) }

pbi_pca.variances = function (pca) {
  imp=summary(pca)$importance[2,]
  var= summary(pca)$rotation[,]^2
  for (i in 1:ncol(var)) { var[,i]=var[,i]*imp[i] }
  return(var)
}

#' \name{pbi_pca.varplot}
#' \alias{pbi$pca.varplot}
#' \alias{pbi_pca.varplot}
#' \title{PCA variance plot to total variances for each component and variable.}
#' \description{
#'     The function provides a PCA matrix plot to show associations 
#'   between PCs and variables. Sown are the squared values, but retaining the 
#'   original sign of the the variances. So the abolute sum of all values should be one.
#' }
#' \usage{pbi_pca.varplot(pca,pcs=10,main="Variance plot",cex.lab=1.5,cex.sym=8, cex.var=1, pch=16,...)}
#' \arguments{
#'   \item{pca}{
#'     pca object which was created using the function _prcomp_.
#'   }
#'   \item{pcs}{
#'     number of the first PC's or vector of PC names to be plotted, default: 10 (or max number of PC's)
#'   }
#'   \item{main}{
#'     title of the plot, default: 'Variance plot'
#'   }
#'   \item{cex.lab}{
#'     character expansion for column and rownames, default: 1.5
#'   }
#'   \item{cex.sym}{
#'     character expansion for the plotting symbols, default: 8
#'   }
#'   \item{cex.var}{
#'     character expansion for the variance values, default: 1
#'   }
#'   \item{pch}{
#'     plotting character for the variances, default: 16 (filled circle)
#'   }
#'   \item{...}{
#'     remaining arguments are delegated to the standard plot function
#'   }
#' }
#' \examples{
#' data(USArrests)
#' pbi_pca.varplot(prcomp(USArrests,scale=TRUE),cex.sym=8)
#' data(swiss)
#' pca=prcomp(swiss,scale=TRUE)
#' round(pbi_pca.variances(pca),3)
#' pbi_pca.varplot(pca,cex.sym=5,cex.var=0.7,cex.lab=1)
#' }
#' \seealso{  See also [pbi](#pbi), [pbi$pca.variances](#pca.variances) }

pbi_pca.varplot = function (pca,pcs=10,main="Variance plot",
                              cex.lab=1.5,cex.sym=8, cex.var=1, pch=16,
                              ...) {
  vars=pbi_pca.variances(pca)
  pcar=summary(pca)$rotation
  if (class(pcs)=="numeric") {
    if (pcs>ncol(vars)) {
      pcs=ncol(vars)
    }
    vars=vars[,1:pcs]
  } else {
    vars=vars[,pcs]
  }
  mt=vars
  yend=nrow(mt)+1
  xend=ncol(mt)+1
  plot(1,type="n",xlab="",ylab="",axes=FALSE,
       xlim=c(-0.5,xend),ylim=c(nrow(mt)+0.35,0),...)
  text(1:(ncol(mt)),0.25,colnames(mt),cex=cex.lab)
  text(-0.5,1:nrow(mt),rownames(mt),cex=cex.lab,pos=4)
  cols=paste("#DD3333",rev(c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD")),sep="")
  cols=c(cols,paste("#3333DD",c(15,30, 45, 60, 75, 90, "AA","BB","CC","DD"),sep=""))
  breaks=seq(-1,1,by=0.1)                  
  mt[pcar<0]=mt[pcar<0]*-1                      
  for (i in 1:nrow(mt)) {
    for (j in 1:ncol(mt)) {
      coli=cut(mt[j,i],breaks=breaks,labels=1:20)
      points(i,j,pch=pch,cex=cex.sym,col=cols[coli])
      text(x=j,y=i,label=sprintf("%0.2f",mt[i,j],2),cex=cex.var)
    }
  }
}

#' \name{pbi_pca.toData}
#' \alias{pbi$pca.toData}
#' \alias{pbi_pca.toData}
#' \title{Transform prcomp pca object  back to data.}
#' \description{
#'   The method allows you transform PCA data back to original data. 
#'   This can be as well used to eliminate some components and then create
#'   data by removing the effect of these components.
#' }
#' \usage{pbi_pca.toData(pca)}
#' \arguments{
#'   \item{pca}{
#'     pca object of class prcomp.
#'   }
#' }
#' \value{return data matrix with the original data.}
#' \examples{
#'   pca=prcomp(iris[,1:4])
#'   head(iris[,1:4])
#'   head(pbi_pca.toData(pca))
#'   # remove effect of first component
#'   pca$rotation[,1]=0
#'   head(pbi_pca.toData(pca))
#' }

pbi_pca.toData = function (pca) {
  if (!is.logical(pca$scale)) {
    data=t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)   
  } else {
    data=t(t(pca$x %*% t(pca$rotation)) + pca$center)
  }
  return(data)
}

#' \name{pbi_pcor}
#' \alias{pbi$pcor}
#' \alias{pbi_pcor}
#' \title{Partial correlation between two variables.}
#' \description{
#'     Calculate partial correlation coefficient of either parametric ("Pearson") 
#'   or non-parametric ("Spearman") statistics corrected for one or more other variables.
#' }
#' \usage{pbi_pcor(x,y,z,method='pearson')}
#' \arguments{
#'   \item{x}{
#'     numeric vector, missing values are allowed
#'   }
#'   \item{y}{
#'     numeric vector, missing values are allowed
#'   }
#'   \item{z}{
#'     numeric vector, matrix or data frame,  missing values are allowed
#'   }
#'   \item{method}{
#'     character string indicating which partial correlation coefficient is to be computed, either "pearson" (default), or "spearman"
#'   }
#' }
#' \value{return partial correlation coefficient between x and y given z.}
#' \examples{
#'   y.data <- data.frame(
#'     hl=c(7,15,19,15,21,22,57,15,20,18),
#'     disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
#'     deg=c(9,2,3,4,1,3,1,3,6,1),
#'      BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00,
#'           4.48e-03,2.10e-06,0.00e+00)
#'   )
#'   # partial correlation between "hl" and "disp" given "deg" and "BC"
#'   pbi_pcor(y.data$hl,y.data$disp,y.data[,c("deg","BC")])
#' }
#' \seealso{  See also [pbi_pcor.test](#pcor.test) }

pbi_pcor = function (x,y,z,method='pearson') {
  r=pbi_pcor.test(x,y,z,method=method)$estimate
  return(r)
}

#' \name{pbi_pcor.test}
#' \alias{pbi$pcor.test}
#' \alias{pbi_pcor.test}
#' \title{Partial correlation test for two variables.}
#' \description{
#'     Calculate partial correlation coefficient and  parametric 
#'   ("Pearson") or non-parametric ("Spearman") 
#'   test statistics for two variables corrected 
#'   for one or more other variables.
#' }
#' \usage{pbi_pcor.test(x,y,z,method='pearson')}
#' \arguments{
#'   \item{x}{
#'     numeric vector, missing values are allowed
#'   }
#'   \item{y}{
#'     numeric vector, missing values are allowed
#'   }
#'   \item{z}{
#'     numeric vector, matrix or data frame,  missing values are allowed
#'   }
#'   \item{method}{
#'     character string indicating which partial correlation coefficient is to be computed, either "pearson" (default), or "spearman"
#'   }
#' }
#' \value{return list with the following components: 
#' 
#' > - _estimate_ - gives the partial correlation coefficient between x and y given z
#'   - _p.value_ - gives the p-value of the test
#'   - _statistics_ - gives the value of the test statistics
#'   - _n_ - gives the number of samples after deleting all the missing samples
#'   - _gn_ - gives the number of given variables
#'   - _method_ - gives the correlation method used}
#'   
#' \examples{
#'   y.data = data.frame(
#'    hl=c(7,15,19,15,21,22,57,15,20,18),
#'    disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
#'    deg=c(9,2,3,4,1,3,1,3,6,1),
#'     BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00,
#'          4.48e-03,2.10e-06,0.00e+00)
#'   )
#'   # partial correlation between "hl" and "disp" given "deg" and "BC"
#'   pbi_pcor.test(y.data$hl,y.data$disp,y.data[,c("deg","BC")])
#' }
#' \seealso{  See also [pcor](#pcor) }

pbi_pcor.test = function (x,y,z,method='pearson') {
  if (is.data.frame(z)) {
    z=as.matrix(z)
  }
  if (is.matrix(z)) {
    df=data.frame(x=x,y=y) 
    for (col in 1:ncol(z)) {
      df=cbind(df,z=z[,col])
      colnames(df)[ncol(df)]=colnames(z)[col]
    }
  } else {
    df=data.frame(x=x,y=y,z=z)
  }
  frmx=formula(paste("x~",paste(colnames(df)[3:ncol(df)],collapse="+"),sep=""))
  frmy=formula(paste("y~",paste(colnames(df)[3:ncol(df)],collapse="+"),sep=""))
  df=na.omit(df)
  if (method=='spearman') {
    df=as.data.frame(apply(df,2,rank))
  }
  #xres=residuals(lm(df[,1]~df[,3]))
  #yres=residuals(lm(df[,2]~df[,3]))
  xres=residuals(lm(frmx,data=df))
  yres=residuals(lm(frmy,data=df))
  
  pr=cor.test(xres,yres,use="complete.obs")
  gn=ncol(df)-2 # number of z
  n=nrow(df)
  statistic <- pr$estimate*sqrt((n-2-gn)/(1-pr$estimate^2))
  p.value <- 2*pnorm(-abs(statistic))
  return(list(estimate=pr$estimate,p.value=pr$p.value,
              conf.int=pr$conf.int,statistic=pr$statistic,
              df=pr$parameter,method=method))
}

#' \name{pbi_prosite2regex}
#' \alias{pbi$prosite2regex}
#' \alias{pbi_prosite2regex}
#' \title{Convert PROSITE patterns to regular expressions.}
#' \description{
#'     This little function converts a PROSITE pattern to a normal regular 
#'   expression. Please note, that for searching in FASTA files those 
#'   regular expressions did work correctly if the sequence pattern
#'   in the FASTA file is splitted over two files, the sequences should be
#'   before fused to single line sequences. The function 
#'   [pbi_searchFasta](#searchFasta) does this.
#' }
#' \usage{pbi_prosite2regex(pattern)}
#' \arguments{
#'   \item{pattern}{
#'     a PROSITE pattern.
#'   }
#' }
#' \value{return a "normal" regular expression.}
#' \examples{
#'   pbi_prosite2regex("A-T-x(0,3)-{ALV}-A")
#' }

pbi_prosite2regex = function (pattern) {
  pattern=gsub("^<","^",pattern)
  pattern=gsub(">$","$",pattern)
  pattern=gsub("\\{([^\\}]+)\\}","[^\\1]",pattern)
  pattern=gsub("\\((.+?)\\)","{\\1}",pattern)
  pattern=gsub("-","",pattern)
  pattern=gsub("x",".",pattern)
  return(pattern)
}

#' \name{pbi_protscale}
#' \alias{pbi$protscale}
#' \alias{pbi_protscale}
#' \title{Calculate protscale moving averages and optionally plot them.}
#' \description{
#'   This functions takes a given input sequence and calculates the hydophobicity averages
#'   over a slising window of length 9 using the Kyte-Doolitle scale. The averaged scores
#'   can be as well plotted.
#' }
#' \usage{pbi_protscale(sequence,plot=FALSE,col='orange')}
#' \arguments{
#'   \item{sequence}{
#'     either a sequence in FASTA format or a simple sequence or text string.
#'   }
#'   \item{plot}{
#'     should the biophysical properties be plotted, default: FALSE
#'   }
#'   \item{col}{
#'     color for the plot 
#'   }
#' }
#' \value{return vector of window averages (invisible).}
#' \examples{
#' fasta="
#' >sp|P04156.1|PRIO_HUMAN
#' MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHG
#' GGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIH
#' FGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVV
#' EQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
#' pbi_protscale(fasta,plot=TRUE,col="light blue")
#' }

pbi_protscale = function (sequence,plot=FALSE,col='orange') {
  # the scores
  kydo=list(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,
            E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,
            F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
  # emtpy result vars
  res=c()
  vals=c()
  main=''
  # split string into lines
  seqs=strsplit(sequence,"\n",perl=TRUE)[[1]]
  for (line in seqs) {
    if (grepl('>',line)) { main=line }
    if (grepl("^\\s*[A-Z]+",line)) {
      # split seq lines into letter
      aas=strsplit(gsub(" ","",line),'')[[1]]
      for (aa in aas) {
        if (aa %in% names(kydo)) {
          vals=c(vals,kydo[[aa]])
        } else {
          stop(paste('Error: unknown AA',aa,'!',sep=''))
        }
      } 
    }
  }
  # now compute the sliding average with window of 9
  for (i in 5:(length(vals)-4)) {
    res=c(res,mean(vals[(i-4):(i+4)]))
  }
  if (plot) {
    plot(5:(length(res)+4),res,type='h',col=col,
         lwd=3, main=main,
         ylab='Score Kyte Doolitte',xlab='Position')
    abline(h=seq(-2,3,by=1),lty=2)
    abline(v=seq(50,250,by=50),lty=2)
  }
  invisible(res)
}

#' \name{pbi_readFasta}
#' \alias{pbi$readFasta}
#' \alias{pbi_readFasta}
#' \title{Read in a (small) fasta file into a list.}
#' \description{
#'     This function can be used to read in a small FASTA file into a listg where the keys are the ids and the values are
#'   are the sequences without line breaks. This code is very slow for large files.
#' }
#' \usage{pbi_readFasta(filename)}
#' \arguments{
#'   \item{filename}{
#'     a sequence file in FASTA format
#'   }
#' }
#' \value{return a list with the sequence ids as keys and the sequences as values.}
#' \examples{
#' fout = file('minions.fasta','w')
#' cat(">Minion1\nSTSTTS\n>Minion2\nTTTTTT\n>Minion3\nSTSTTT\n",file=fout)
#' cat(">Minion4\nSTTTTT\n>Minion5\nSSTTTT\n>Minion6\nSSSTST\n",file=fout)
#' cat("\n>Minion7\nSSSSTT\n", file=fout)
#' close(fout)
#' seq=pbi_readFasta("minions.fasta")
#' seq[7]
#' M=adist(unlist(seq))
#' M
#' plot(hclust(as.dist(M)))
#' }
#' \seealso{  See also [pbi](#home), [pbi_searchFasta](#searchFasta) }

pbi_readFasta <- function (filename) {
  res=list()
  if (!file.exists(filename)) {
    stop(paste("Error: file",filename,"does not exists!"))
  }
  fin  = file(filename, "r")
  seq=""
  while(length((line = readLines(fin,n=1)))>0) {
    if (grepl("^>([^ ]+)",line)) {
      if (seq != "") {
        res[id]=seq
      }
      id=gsub("^>([^  ]+).*$","\\1",line)
      seq=""
    } else {
      seq=paste(seq,line,sep="")
    }
  }
  close(fin)
  res[id]=seq
  return(res)
}

#' \name{pbi_readPepinfo}
#' \alias{pbi$readPepinfo}
#' \alias{pbi_readPepinfo}
#' \title{Read data from the EMBOSS pepinfo tool.}
#' \description{
#'   This is a function to visualize the biophysical properties of a protein
#'   using the output of the EMBOSS pepinfo tool which can be accessed online
#'   at [https://www.ebi.ac.uk/Tools/seqstats/emboss_pepinfo/](https://www.ebi.ac.uk/Tools/seqstats/emboss_pepinfo/). After submitting your sequence you have to use the file at "Result Files->Tool Output" the file ending with ".output".
#' }
#' \usage{pbi_readPepinfo(file,region=NULL)}
#' \arguments{
#'   \item{file}{
#'     the result file from the EMBOSS pepinfo file.
#'   }
#'   \item{region}{
#'     string matching a certain region such as "Doolittle", if not given all available regions will be shown,  default: NULL
#'   }
#' }
#' \value{return dataframe with the columns: Position, Aminoacid, Result.}
#' \examples{
#' pepfile=system.file("files/pepinfo-spike-sars2.txt",package="pbi")
#' if (file.exists(pepfile)) {
#'   pbi_readPepinfo(pepfile)
#'   res=pbi_readPepinfo(pepfile,region="Doolittle")
#'   plot(res$Result ~ res$Position,type="h",col="orange")
#' }
#' }

pbi_readPepinfo = function (file,region=NULL) {
  data=read.table(file,sep="\n",header=FALSE,
                  blank.lines.skip=TRUE)
  idx=grep("Printing out",data$V1)
  res=data.frame()
  if (is.null(region)) {
    print(data[idx,])
    return()
  }
  start=intersect(grep(region,data$V1),idx)
  print(start)
  j=0
  for (i in (start+1):nrow(data)) {
    if (length(grep("Printing out",data[i,'V1']))>0) {
      break
    } else if (length(grep("Created emboss",data[i,'V1']))>0) {
      break
    }else if (length(grep("Position",data[i,'V1']))>0) {
      next
    } else {
      j=j+1
      row=strsplit(gsub(" +"," ",data[i,'V1'])," ")[[1]]
      #print(row)
      if (j==1) {
        res=data.frame(Position=as.integer(row[2]),
                       Aminoacid=row[3],
                       Result=as.numeric(row[4]))
      } else {
        res=rbind(res,data.frame(Position=as.integer(row[2]),
                                 Aminoacid=row[3],
                                 Result=as.numeric(row[4])))
      }
    }
  }
  return(res)
}

#' \name{pbi_report.chisq.test}
#' \alias{pbi$report.chisq.test}
#' \alias{pbi_report.chisq.test}
#' \title{Return a formatted text string for reporting a chisq.test.}
#' \description{
#'     Return a formatted text string for reporting a chisq.test.
#' }
#' \usage{pbi_report.chisq.test(tab)}
#' \arguments{
#'   \item{tab}{
#'     a contigency table.
#'   }
#' }
#' \value{return formatted text string for reporting a chisq.test in a LaTeX/Sweave document.}
#' \examples{
#'   azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2))
#'   rownames(azt)=c("AZT","Placebo")
#'   colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress")
#'   pbi_report.chisq.test(azt)
#' }
#' \seealso{  See also [pbi](#home), [pbi_report.pval](#report.pval) }

pbi_report.chisq.test <- function (tab) {
  ct=prop.test(tab)
  return(paste('$\\\\chi^2$(',ct$parameter[[1]],',',
               sum(tab),')=',round(ct$statistic,2),
               ', p',pbi_report.pval(ct$p.value), 
               ', Cohens $w$=',round(pbi_cohensW(tab),2),sep=''))
}

#' \name{pbi_report.conf.int}
#' \alias{pbi$report.conf.int}
#' \alias{pbi_report.conf.int}
#' \title{Return a formatted text string for reporting a confidence interval. }
#' \description{
#'     Return a formatted text string for reporting a confidence interval. 
#' }
#' \usage{pbi_report.conf.int(ci,round=2)}
#' \arguments{
#'   \item{ci}{
#'     a confidence interval consisting of two numerical values
#'   }
#'   \item{round}{
#'     rounding digits, default: 2
#'   }
#' }
#' \value{return a formatted text string for reporting a confidence interval to be included in a LaTeX/Sweave document.}
#' \examples{
#'   azt=as.table(matrix(c(76,399,129,332), byrow=TRUE,ncol=2))
#'   rownames(azt)=c("AZT","Placebo")
#'   colnames(azt)=c("DiseaseProgress", "NoDiseaseProgress")
#'   pbi_report.conf.int(prop.test(azt)$conf.int)
#' }
#' \seealso{  See also [pbi](#home), [pbi_report.pval](#report.pval) }

pbi_report.conf.int <- function (ci,round=2) {
  return(paste('[',round(ci[1],round),',',
               round(ci[2],round),']',sep=''))
}

#' \name{pbi_report.pval}
#' \alias{pbi$report.pval}
#' \alias{pbi_report.pval}
#' \title{Return a p-value of reporting, either giving the three alpha thresholds, 
#'   <0.05, <0.01, or <0.001 or using the star syntax. }
#' \description{
#'     Return a p-value of reporting, either giving the three alpha thresholds, 
#'   <0.05, <0.01, or <0.001 or using the star syntax. 
#' }
#' \usage{pbi_report.pval(p.val,star=FALSE)}
#' \arguments{
#'   \item{p.val}{
#'     a numerical p-value.
#'   }
#'   \item{star}{
#'     boolean, should the one-three star syntax be used, default: FALSE.
#'   }
#' }
#' \value{return pvalue in shown in alpha-threshold or star syntax. If the p-value is not significant, either the value or any empty string is returned if the star syntax is used.}
#' \examples{
#'   report.pval = pbi_report.pval
#'   report.pval(1/10000)
#'   report.pval(1/10000,star=TRUE)
#'   report.pval(0.02,star=TRUE)
#'   report.pval(0.12,star=TRUE)
#'   report.pval(c(0.001,0.01,0.3,0.02))
#' }
#' \seealso{  See also [pbi](#home), [pbi_report.conf.int](#report.conf.int) }

pbi_report.pval <- function (p.val,star=FALSE) {
  if (length(p.val) > 1) {
    return(as.character(lapply(p.val,pbi_report.pval)))
  }
  if (p.val <0.001 & star) {
    return('***')
  } else if (p.val <0.001) {
    return('<0.001')
  } else if (p.val <0.01 & star) {
    return('**')
  } else if (p.val <0.01) {
    return('<0.01')
  } else if (p.val <0.05 & star) {
    return('*')
  } else if (p.val <0.05) {
    return('<0.05')
  } else if (star) {
    return("")
  } else {
    return(sprintf("%.2f",p.val))
  }   
}   

#' \name{pbi_sem}
#' \alias{pbi$sem}
#' \alias{pbi_sem}
#' \title{Calculate the standard error of the mean for a given numerical vector.}
#' \description{
#'     Calculate the standard error of the mean for a given numerical vector.
#' }
#' \usage{pbi_sem(x,na.rm=FALSE)}
#' \arguments{
#'   \item{x}{
#'     a numerical vector.
#'   }
#'   \item{na.rm}{
#'     logical vecor indicating if missing values should be removed, default: FALSE
#'   }
#' }
#' \value{return computed standard error of the mean.}
#' \examples{
#'   pbi_sem(rnorm(50,mean=10,sd=3))
#'   pbi_sem(rnorm(1000,mean=10,sd=3))
#' }
#' \seealso{  See also [pbi](#home), [pbi_cv](#cv) }

pbi_sem <- function(x,na.rm=FALSE) {
  sd(x,na.rm=na.rm)/sqrt(length(x[!is.na(x)])) 
}

#' \name{pbi_searchFasta}
#' \alias{pbi$searchFasta}
#' \alias{pbi_searchFasta}
#' \title{Search a FASTA file with a regular expression.}
#' \description{
#'     This function searches FASTA files by removing line breaks within the
#'   sequence belonging to the same ID.
#' }
#' \usage{pbi_searchFasta(filename,pattern)}
#' \arguments{
#'   \item{filename}{
#'     a sequence file in FASTA format
#'   }
#'   \item{pattern}{
#'     a standard regular expression which will be applied on the sequence for the ID.
#'   }
#' }
#' \value{return the matching ID's.}
#' \examples{
#'   # the pattern splits over two lines for the first sequence
#'   pbi_searchFasta(filename=system.file("files/human-tRNAs.fasta",package="pbi"),
#'      pattern="GACGC{5}AT{2}CTCT")
#' }
#' \seealso{  See also [pbi](#home), [pbi_tkregex](#tkregex), [pbi_prosite2regex](#prosite2regex) }

pbi_searchFasta <- function(filename,pattern) {
  if (!file.exists(filename)) {
    stop(paste("Error: file",filename,"does not exists!"))
  }
  fin  = file(filename, "r")
  ids=c()
  seq=""
  while(length((line = readLines(fin,n=1)))>0) {
    if (grepl("^>([^ ]+)",line)) {
      if (seq != "") {
        if (grepl(pattern,seq)) {
          ids=c(ids,id)
        }
      }
      id=gsub("^>([^  ]+).*$","\\1",line)
      seq=""
    } else {
      seq=paste(seq,line,sep="")
    }
  }
  close(fin)
  if (grepl(pattern,seq)) {
    ids=c(ids,id)
  }
  return(ids)
}

#' \name{pbi_text2fasta}
#' \alias{pbi$text2fasta}
#' \alias{pbi_text2fasta}
#' \title{Convert any Text file to a FASTA file using only Aminoacid letters.}
#' \description{
#'     This function loops over the given directory and converts text files to
#'   a FASTA file format. Only uppercase letters are used for the output file.
#' }
#' \usage{pbi_text2fasta(dir,pattern="*",outfile="stdout")}
#' \arguments{
#'   \item{dir}{
#'     a directory
#'   }
#'   \item{pattern}{
#'     regular expression for the files which should be analysed
#'   }
#'   \item{outfile}{
#'     the output filename, default: 'stdout'
#'   }
#' }
#' \value{return the number of sequences.}
#' \examples{
#'   pbi_text2fasta(dir=".",pattern='*.R$',outfile="test.fasta")
#'   pbi_file.head("test.fasta",n=3)
#' }
#' \seealso{  See also [pbi](#home), [pbi_searchFasta](#searchFasta) }

pbi_text2fasta <- function (dir,pattern="*",outfile="stdout") {
  if (outfile == "stdout") {
    out = stdout()
  } else {
    out = file(outfile,"w")
  }
  for (file in list.files(dir,full.names=TRUE,
                          pattern=pattern)) {
    fin = file(file,'r')
    cat('>',gsub(".+/","",file),'\n',sep="",file=out)
    seq = ""
    while(length((
      line = readLines(fin,n=1,warn=FALSE)))>0) {
      line = gsub("#+","#",line)
      line = gsub("[0-9]","X",line)
      line = gsub("[\\s\\t]","Y",line)
      line = gsub("[^A-Za-z]","Z",line)            
      seq=paste(seq,line,sep="")
    }
    seq = toupper(seq)
    close(fin)
    cat(gsub("([A-Z]{70})","\\1\n",seq),"\n",sep="",file=out)
  }
  close(out)
}

#' \name{pbi_tkregex}
#' \alias{pbi$tkregex}
#' \alias{pbi_tkregex}
#' \title{Graphical user interface for testing regular expressions.}
#' \description{
#'   This function provides a graphical user interface to execute regular
#'   expressions either as PROSITE patterns or a standard regular expressions
#'   on the text entered or pasted into the text area. Just copy and paste 
#'   your text into the text area below of the two entry fields. 
#'   In the entry fields you can write your regular expression either as normal regular
#'   expression in the bottom field or as PROSITE pattern in the top entry.
#'   Thereafter press ENTER, your found patterns will be highlighted in red and 
#'   the number of found patterns will be shown as well.
#' }
#' \usage{pbi_tkregex()}
#' \examples{
#'   \dontrun{ 
#'     pbi_tkregex()
#'   }
#' }
#' \seealso{  See also [pbi](#home), [pbi_prosite2regex](#prosite2regex) }

pbi_tkregex = function () {
  source(system.file("scripts/tkregex.R",package="pbi"))
}

#' \name{pbi_wordFreq}
#' \alias{pbi$wordFreq}
#' \alias{pbi_wordFreq}
#' \title{Count the number of words with a given length in a sequence.}
#' \description{
#'   The function creates a sliding window of length `wlength` over the given 
#'   text or sequence string creating words of length `wlength`. 
#'   The number how often a certain word appears is counted. 
#'   So - `AABAAC` - contains the words: `AA`, `AB`, `BA`, `AA` and `AC`.
#' }
#' \usage{pbi_wordFreq(seq,wlength=2)}
#' \arguments{
#'   \item{seq}{
#'     sequence as text string or vector of text strings
#'   }
#'   \item{wlength}{
#'     word size, default: 2
#'   }
#' }
#' \value{return list object with words as keys and counts as values or data frame in case input is a vector.}
#' \examples{
#'   seq="AAABBBCCCDEFAABBCCDD"
#'   unlist(pbi_wordFreq(seq))
#'   unlist(pbi_wordFreq(seq,wlength=4))
#'   minions=read.table(text='
#'   Minion1      STSTTS
#'   Minion2      TTTTTT
#'   Minion3      STSTTT
#'   Minion4      STTTTT
#'   Minion5      SSTTTT
#'   Minion6      SSSTST
#'   Minion7      SSSSTT
#' ',row.names=1)
#'   min=pbi_wordFreq(minions[,1],wlength=2)
#'   rownames(min)=rownames(minions)
#'   min
#' }

pbi_wordFreq = function (seq,wlength=2) {
  if (length(seq) > 1) {
    df=NULL
    for (i in 1:length(seq)) {
      res=pbi_wordFreq(as.character(seq[i]),wlength=wlength)
      ul=unlist(res)
      if (class(df[1])=="NULL") {
        df= t(data.frame(ul))
      } else {
        nnames= setdiff(names(ul),colnames(df))
        if (length(nnames) > 0) {
          for (n in nnames) {
            df=cbind(df,ncol=rep(0,nrow(df)))
            colnames(df)[ncol(df)]=n
          }
        }
        df=rbind(df,rep(0,ncol(df)))
        df[nrow(df),names(ul)]=as.vector(ul)
      }
      
      
    }
    rownames(df)=1:nrow(df)
    return(df)
  } else {
    res=list()
    for (i in 1:(nchar(seq)-wlength+1)) {
      s=substr(seq,i,i+wlength-1)
      if (!is.null(res[[s]])) {
        res[[s]]=res[[s]]+1
      } else {
        res[[s]]=1
      }
    }   
    return(res)
  }
}

#' \name{pbi_wininstall}
#' \alias{pbi$wininstall}
#' \alias{pbi_wininstall}
#' \title{Make an executable Batch script on Windows for R applications within the users PATH.}
#' \description{
#'   The function two files in the users PATH for executables on Windows, a BATCH file and a script file
#'   both in the same folder and with the same file prefix. This allows the user to directly execute a Rscript 
#'   by pressing Win-R as a shortcut and the then entering the application name assume you have a file: _hw.R_
#' }
#' \usage{pbi_wininstall(filename="")}
#' \arguments{
#'   \item{filename}{
#'     optional filename which should be installed directly without asking for the filename, default: ""
#'   }
#' }
#' \examples{
#'   \dontrun{
#'   pbi_wininstall("hw.R")
#' }
#' }

pbi_wininstall = function (filename="") {
  getUserPath = function () {
    p=""
    for (path in sort(unique(strsplit(Sys.getenv("PATH"),";")[[1]]))) {
      if (grepl("Users.+WindowsApps",path)) {
        return(path)
      } else if (grepl("C:.*Users",path)) {
        p=path
      }
    }
    return(p)
  }
  
  drive=substr(R.home(),1,1)
  batchcode="@echo off
dir \"DRIVE:\\Program Files\\Rscript.exe\" /s /b |  findstr x64 >\"%~dp0\\Rversion.txt\"
set /p RS=<\"%~dp0\\Rversion.txt\"
\"%RS%\" \"%~dpn0.R\" %1 %2 %3 %4 %5 %6 %7 %8 %9
"  
  rscriptcode="@echo off
dir \"DRIVE:\\Program Files\\Rscript.exe\" /s /b |  findstr x64 >\"%~dp0\\Rversion.txt\"
set /p RS=<\"%~dp0\\Rversion.txt\"
\"%RS%\" %1 %2 %3 %4 %5 %6 %7 %8 %9
"  
  batchcode=gsub("DRIVE",drive,batchcode)
  rscriptcode=gsub("DRIVE",drive,rscriptcode)
  
  if (filename == ""){
    filename=tclvalue(tcltk::tkgetOpenFile(filetypes="{{R files} {*.r *.R}} {{All files} {*.*}}"))
  }
  path=getUserPath()
  if (filename != "") {
    #print(basename(filename))
    bname=paste(gsub(".[Rr]$","",basename(filename)),".bat",sep="")
    #print(bname)
    outfile=file(file.path(path,bname),"w")
    cat(batchcode,file=outfile)
    close(outfile)
    file.copy(filename,path,overwrite = TRUE)
    #print(file.exists(file.path(path,bname)))
    print(paste("Install of ",filename,"was successful installed!"))
    print(paste("You can now call the executable with:",bname))
    print(path)
    Sys.sleep(5)
  }
  rscript=file.path(path,"Rscript.bat")
  if (!file.exists(rscript)) {
    rsout=file(rscript,"w")
    cat(rscriptcode,file=rsout)
    close(rsout)
    print(paste("Rscript.bat created at:",path))
    Sys.sleep(5)
    
  }
}

#' \name{pbi_xyplot}
#' \alias{pbi$xyplot}
#' \alias{pbi_xyplot}
#' \title{Improved xyplot which as well displays a grid and the correlation coefficient.}
#' \description{
#'   This is just a simple illustrative function to demonstrate how a standard
#'   R function can be modified and extended using different type of arguments
#'   overwriting default settings and using the ellipsis and delegating 
#'   remaining arguments to the default function
#' }
#' \usage{pbi_xyplot(x,y,col="blue",pch=19,grid=TRUE,xlab=NULL,ylab=NULL,ellipse=FALSE,ell.fill=FALSE,show.r=TRUE,...)}
#' \arguments{
#'   \item{x}{
#'     vector with  numerical values, or a matrix or data frame
#'   }
#'   \item{y}{
#'     vector with  numerical values, can be NULL if x is matrix or data frame, default: NULL
#'   }
#'   \item{col}{
#'     color for the plotting symbols
#'   }
#'   \item{pch}{
#'     plotting character, default: 19 (filled circle)
#'   }
#'   \item{grid}{
#'     should be a grid drawn, default: TRUE
#'   }
#'   \item{xlab}{
#'     xlabel, default: NULL (x)
#'   }
#'   \item{ylab}{
#'     ylabel, default: NULL (y)
#'   }
#'   \item{ellipse}{
#'     draw an ellipse for 85% and 95% confidence intervals using the given colors, default: FALSE
#'   }
#'   \item{ell.fill}{
#'     should the ellipse be filled for the 85% confidence intervall,  default: FALSE
#'   }
#'   \item{show.r}{
#'     should the Pearson correlation been shown on top, default: TRUE
#'   }
#'   \item{...}{
#'     remaining arguments delegated to the standard plot function
#'   }
#' }
#' \examples{
#' data(iris)
#' par(mfrow=c(2,2),mai=c(0.8,0.8,0.8,0.1))
#' pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length,
#'     xlab="Sepal.Width",ylab="Sepal.Height",
#'     col=as.numeric(iris$Species)+1)
#' pbi_xyplot(iris$Petal.Width,iris$Petal.Length,
#'     xlab="Petal.Width",ylab="Petal.Height",
#'     col=as.numeric(iris$Species)+1)
#' legend("bottomright",fill=2:4,legend=levels(iris$Species))
#' pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length,
#'     xlab="Sepal.Width",ylab="Sepal.Height",
#'     col=as.numeric(iris$Species)+1,ellipse=TRUE)
#' pbi_xyplot(iris$Sepal.Width,iris$Sepal.Length,
#'     xlab="Sepal.Width",ylab="Sepal.Height",
#'     col=as.numeric(iris$Species)+1,ellipse=TRUE,ell.fill=TRUE)
#' }

pbi_xyplot <- function(x,y,col="blue",pch=19,grid=TRUE,
                         xlab=NULL,ylab=NULL,ellipse=FALSE,
                         ell.fill=FALSE,show.r=TRUE,...) {
  if (is.matrix(x) | is.data.frame(x)) {
    y=x[,2]
    x=x[,1]
  }
  if (is.null(xlab)) {
    xlab="x"
  } 
  if (is.null(ylab)) {
    ylab="y"
  } 
  plot(x,y,col=col,pch=pch,xlab=xlab,ylab=ylab,...)
  if (grid) {
    grid (NULL,NULL, lty = 3, col = "grey30")
  }
  if (ellipse) {
    if (length(col) == 1) {
      col=rep(col,length(x))
    }
    if (length(col)!= length(x)) {
      stop("colors must have same length as data points")
    }
    ell.col=col
    i=1
    for (cl in names(table(ell.col))) {
      C=cov(data.frame(x=x[ell.col==cl],y=y[ell.col==cl]))    # Covarianz-Matrix C bestimmen
      d85=qchisq(0.85, df = 2)     # 85% - Faktor , um die Ellipse zu skalieren
      M=colMeans(data.frame(x=x[ell.col==cl],y=y[ell.col==cl])) #   Mittelwerte (Zentrum) des Clusters
      el=cluster::ellipsoidPoints(C, d85, loc=M)  # Ellipsen-Punkte aus C und M berechnen
      if (ell.fill) {
        colfill=paste(rgb(t(col2rgb(cl))/255),"33",sep="")
        polygon(el,col=colfill,border=NA)
        i=i+1
        next
      }
      lines(el,col=cl,lwd=1.5,lty=2)    #  Ellipse als geschlossene Linies zeichnen
      d95=qchisq(0.95, df = 2)     # 85% - Faktor , um die Ellipse zu skalieren
      el=cluster::ellipsoidPoints(C, d95, loc=M)  # Ellipsen-Punkte aus C und M berechnen
      lines(el,col=cl,lwd=1.5,lty=1)    #  Ellipse als geschlossene Linies zeichnen                        
    }
  }
  if (show.r) {
    ct=cor.test(x,y)
    mtext(paste("r =",round(ct$estimate,3),
                pbi_report.pval(ct$p.value,star=TRUE)),
          side=3)
  }
}

#' \name{testprint}
#' \alias{testprint}
#' \title{ print a test message }
#' \description{
#'     A test function which prints a string.
#' }
#' \usage{ testprint(txt) }
#' \arguments{
#'   \item{txt}{ some value to print  }
#' }
#' \examples{
#'     testprint("Hello World!")
#'     # sample file use
#'    dec=read.table(file.path(system.file(package="pbi"),"files","decathlon.tab"))
#'    head(dec)
#' }

testprint <- function (txt) {
  print(txt)
}

## Functions or variables starting with uppercase letters
## will be per default not export, they can be used as 
## internal package functions not accessible by the user
## of the package

Hidden = function (x) {
    return(x+1)
}


###  place startup loads here
.onLoad <- function(libname, pkgname) {
    # to show a startup message
    # example on how to initialize a Tcl package 
    # tcltk::.Tcl(paste("lappend auto_path",file.path(system.file(package="pbi"),"pantcl", "lib")))
    # tcltk::.Tcl("package require tclfilters")
    # tools::vignetteEngine("pantcl",package=pkgname,weave=pantcl,tangle=ptangle)
}
