#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

## Install DropletUtils
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("DropletUtils")
BiocManager::install("BiocParallel")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.8")

library(BiocManager)
BiocManager::valid()

if ( ! library('DropletUtils', character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }
if ( ! library('BiocParallel', character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }


## Install OptParse
install.packages('optparse', dependencies=TRUE)

if ( ! library('optparse', character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }