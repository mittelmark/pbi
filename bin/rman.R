#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)

extractRd <- function (args) {
    if (!file.exists("DESCRIPTION")) {
        stop("Script must be run in the package folder!")
    } 
    pkgname=basename(getwd())
    if (pkgname != "Rpkg" & file.exists("man/Rpkg-package.Rd")) {
        replaceRpkg()
    }
    for (file in args) {
        if (!file.exists(file)) {
            stop(paste("Error: File",file,"does not exists"))
        }
        fin  = file(file, "r")
        flag = FALSE
        while(length((line = readLines(fin,n=1)))>0) {
            if (grepl("^#' *\\\\name\\{.+\\}",line)) {
                      name = gsub("\\$","_",gsub("^#' *\\\\name\\{(.+)\\}.*","\\1",line))
                                  print(paste("name = '",name,"'",sep=''))
                flag = TRUE
                mandir = sub("R/.+","man/",file)
                fout = file(paste(mandir,name,".Rd",sep=""),"w")
                
            } else if (flag & (grepl('^[a-zA-Z" ].*',line) | grepl("^$",line))) {
                flag = FALSE
                close(fout) 
                next
            } 
            if (flag) {
                line = sub("^#' ?","",line)
                cat(line,"\n",file=fout)
            }
        }
        close(fin)
    }
}

replaceRpkg <- function (pkgname=NULL) {
    if (is.null(pkgname)) {
        if (!file.exists("DESCRIPTION")) {
            stop("Script must be run in the package folder!")
        } 
        pkgname=basename(getwd())
    }
    for (file_path in c("tests/add-tests.R","R/add.R","DESCRIPTION")) {
        file_content <- readLines(file_path, warn = FALSE)
        mod = gsub("Rpkg",pkgname,file_content)
        writeLines(mod, file_path)
    }
    file.remove("man/Rpkg-package.Rd")
}
 
main <- function (argv) {
    extractRd(argv)
}
if (sys.nframe() == 0L && !interactive()) {
    main(commandArgs(trailingOnly=TRUE))
}
