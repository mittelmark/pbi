\name{pbi_domainPlot} 
\alias{pbi$domainPlot} 
\alias{pbi_domainPlot} 
\title{Use the mydomain website of prosite to create a protein domain plot.} 
\description{ 
  This function can be used to draw a plot representing protein domains into a standard R plot. 
  The data are send to the website [https://prosite.expasy.org/cgi-bin/prosite/mydomains/](https://prosite.expasy.org/cgi-bin/prosite/mydomains/), a URL is created and the image from this url is downloaded to the local file system. The filename is a CRC32 digest of the URL so allowing to cache as well the downloaded results. 
} 
\usage{pbi_domainPlot(domains,ranges,sites,length=1000,hscale=1,cache=TRUE,plot=TRUE,cex=0.8,...)} 
\arguments{ 
  \item{domains}{ 
    list where names are the domain names and the values are four integers: start, end, shape (1:6), color (1:4) 
  } 
  \item{ranges}{ 
    list with vectors of three integers: start, end, type (0:1) 
  } 
  \item{sites}{ 
    list with vectors of two integers: position, type (0:1) 
  } 
  \item{length}{ 
    length of the protein, default: 1000 
  } 
  \item{hscale}{ 
    horizontal scaling factor (between 0.1 and 2.0), default: 1.0  
  } 
  \item{cache}{ 
    should the image locally cached, default: TRUE 
  } 
  \item{plot}{ 
    should the image plotted on a R plot devices, default: TRUE 
  } 
  \item{cex}{ 
    character expansion for the domain text and sites if plot is TRUE, default: 0.8 
  } 
  \item{...}{ 
    any argument forwarded to the plot function if plot is TRUE 
  } 
} 
\value{return a URL or a local filename in case of cache=TRUE, if plot is TRUE NULL is returned invisible} 
\examples{ 
url=pbi_domainPlot( 
    domains=list(MYDOM1=c(100,200,2,1),MYDOM2=c(300,500,3,2)), 
    ranges=list(a=c(190,500,1)), 
    sites=list(a=c(150,1)),hscale=2.0,plot=FALSE) 
print(url) 
} 
