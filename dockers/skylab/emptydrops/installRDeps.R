#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

## Install OptParse
install.packages('optparse', dependencies=TRUE, repos="http://R-Forge.R-project.org")

if ( ! library('optparse', character.only=TRUE, logical.return=TRUE) ) {
        quit(status=1, save='no')
    }

## Install DropletUtils
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DropletUtils", version = "3.8")
