#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

## Install OptParse
install.packages('remotes')
remotes::install_version('optparse', '1.6.6')

## Install DropletUtils
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DropletUtils", version = "3.8")
