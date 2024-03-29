<!-- README.md is generated from README.Rmd. Please edit that file -->


# gwhap

<!-- badges: start -->
<!-- badges: end -->

gwhap, the Genome-Wide HAPlotype analysis package

The gwhap software is related to the following paper : https://www.nature.com/articles/s41431-021-00827-8 

## Installation with a singularity container DURING development cycles

1.  get a singularity that has devetools and IRkernels.
2.  get a clone of the gwhap gitlab project
3.  organize your directories to reflect the following envir varaibles (where is singularity and in which dir is gwhap.):

```bash
SINGU=/home/USERNAME/volatile/singularity/brainomics-ubuntu.simg
MY_HOME=/home/USERNAME/volatile/PROJECT
```

Then run :
```bash
singularity shell --home $MY_HOME:/home/jovyan $SINGU
```

From **within** the singularity
```bash
export R_LIBS=~/R
```

Then launch R
```r
setwd("gwhap")
# this has to be replayed EACH time you modify the code
###############################
detach("package:gwhap", unload=TRUE)

devtools::document(roclets=c('rd', 'collate', 'namespace'))
devtools::build_vignettes()
devtools::build()


install.packages("../gwhap_0.2.tar.gz", repos=NULL)
library(gwhap)

###############################
# now run an exmaple
toy_to_run = system.file("example", "toy.R",package="gwhap")
source(toy_to_run)
```

or 

```bash
R -e "setwd(\"gwhap\");devtools::document(roclets=c('rd', 'collate', 'namespace'))"
R CMD build --no-build-vignettes gwhap ; R CMD INSTALL gwhap_0.2.tar.gz
```


## Installation

If run from within a singularity container:
```bash
export R_LIBS=<PAHTH into the container>/R
```

```R
install.packages("gwhap.xxx.tar.gz",repos=NULL)
```

## Documentation

package documentation is available after installation:
(might need to reload R/rstudio session)
```{r package documentation, echo=TRUE}
package?gwhap 
```
or once it is loaded :
```{r package documentation 2, echo=TRUE}
library(gwhap)
?gwhap 
```
description and full documentation :
```{r package full documentation, echo=TRUE}
library(help="gwhap")

```
