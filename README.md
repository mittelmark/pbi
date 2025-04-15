# pbi - R package for Practical Bioinformatics 

![](../../actions/workflows/r.yml/badge.svg) [![Release](https://img.shields.io/github/v/release/mittelmark/pbi.svg?label=current+release)](https://github.com/mittelmark/pbi/releases) [![license](https://img.shields.io/badge/license-MIT-lightgray.svg)](https://opensource.org/license/MIT)


This is a package  for the  course  Practical  Bioinformatics  for the  Master
course Molecular Biology and Biochemistry at the University of Potsdam.
The package contains functions and data used in the lectures and exercises of this course.

## Installation

To install the package you can either use the package  remotes to install from
the current github sources or you download and install the latest release.

First let's show how to install the latest version directly from Gihub:

```
> install.packages("remotes")
> library(remotes)
> remotes::install_github("https://github.com/mittelmark/pbi")
```

To install instead the latest release version which should be more stable use the following command:

```
> install.packages("https://github.com/mittelmark/pbi/releases/download/v0.1.2/pbi_0.1.2.tar.gz",repos=NULL)
```

## Package documentation

The package documentation for version 0.1.0 can be seen here: [pbi-manual.pdf](https://github.com/mittelmark/pbi/files/14588148/pbi-manual.pdf)

The vignette can be viewed after installing the package by writing into the terminal this: `vignette("pbi-tutorial")`


## Author and Copyright

Author: Detlef Groth, University of Potsdam, Germany

Co-authors: 

- Vera Burdova, University of Potsdam, Germany
- Jillian Hoffmann, University of Potsdam, Germany
- Lisa Banu Kurt, University of Potsdam, Germany

License: MIT License see the file [LICENSE](LICENSE) for details.

## Bug reporting

In case of bugs and suggestions, use the [issues](../../issues) link on top.

