## function documentation is directly embedded within
## the R code
## you can automatically extract the documentation using the bin/rman.R script
## by executing: using this like Rscript bin/rman.R R/add.R
## or from within the R-console:
## > setwd("package-folder")
## > source('bin/rman.R')
## > extractRd(list.files("R",pattern="*.R$",full.names=TRUE))

## package-documentation
#' \name{Rpkg-package}
#' \alias{Rpkg-package}
#' \title{Rpkg package - package template for R packages}
#' \description{The Rpkg package can be used as a template to create new packages from scratch.}
#' \details{Some more details:
#' The following list of objects and/or functions are available:
#' \describe{
#' \item{\link[Rpkg:add]{add(x,y)}}{an illustrative add function}
#' }
#' }
#' \examples{
#' library(Rpkg)
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

#' \name{add}
#' \alias{add}
#' \title{ add two numbers }
#' \description{
#'     A initial starting function which can add two numbers.
#' }
#' \usage{ add(x, y) }
#' \arguments{
#'   \item{x}{
#'     numercial value
#'   }
#'   \item{y}{
#'     numerical value
#'   }
#' }
#' \details{
#'     Some more details ...
#' }
#' \value{return the sum of x and}
#' \examples{
#'     add(2,3)
#'     x=2
#'     y=3
#'     add(x,y)
#' }
#' \seealso{  See also \code{\link{testprint}} }
#' 
#' \keyword{ arith }

add <- function (x, y) {
    x + y
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
#'    dec=read.table(file.path(system.file(package="Rpkg"),"files","decathlon.tab"))
#'    head(dec)
#' }
#' \seealso{  See also \code{\link{add}} }
#' \keyword{ arith }

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
    # tcltk::.Tcl(paste("lappend auto_path",file.path(system.file(package="Rpkg"),"pantcl", "lib")))
    # tcltk::.Tcl("package require tclfilters")
    # tools::vignetteEngine("pantcl",package=pkgname,weave=pantcl,tangle=ptangle)
}
