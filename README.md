# mvmstan
mvmstan simplifies the specification of multivariate mixture models in Stan (mc-stan.org).

# Installation

## 1. Install RGraphviz from bioconductor
 	source("http://bioconductor.org/biocLite.R")
 	biocLite("Rgraphviz")

## 2. Install mvmstan in R
 	install.packages("devtools")
 	library(devtools)
 	install_github("plogacev/mvmstan")
or
 - download the zip archive
 - install the package from the zip file

## 3. Install a Graphviz viewer

### On Linux
 	sudo apt-get install xdot

### On Mac OS X
- A graphviz viwer can be installed from http://www.graphviz.org/Download_macos.php
