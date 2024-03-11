\name{pbi_package.deps} 
\alias{pbi$package.deps} 
\alias{pbi_package.deps} 
\title{Package(s) which are required by the given package name} 
\description{The function helps in identifying the nested package dependencies for a given package.} 
\usage{ pbi_package.deps(pkgName,mode='all') } 
\arguments{ 
  \item{ pkgName }{ 
    an package name given as text string. 
  } 
  \item{ mode }{ 
    which package names to return, the following modes are available: 
    'all' - all required packages,   
    'install' - not yet installed packages,  
    'nonbase' - packages not in the standard R installation 
  } 
} 
\value{return list of required packages.} 
\examples{ 
  \dontrun{ 
    pbi_package.deps('igraph',mode='nonbase') 
    pbi_package.deps('igraph',mode='all') 
 } 
} 
\seealso{  
\code{\link[pbi:pbi-package]{pbi-package}}, \code{\link[pbi:pbi-class]{pbi-class}} 
} 
