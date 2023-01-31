#!/usr/bin/env Rscript

## Description: Wrapper for running emptyDrops on a sparse matrix, provides
##    command line arguments for emptyDrops and I/O handling and
##    error trapping.

## Load optparse
library('optparse')

## Define helper functions

#' Prints a message to stderr and exits R with error code 1
#' @param msg message to standard error
errorExit <- function(msg) {
    cat(msg,file=stderr())
    quit(save='no',status=1)
}

#' Prints a message with cat only if verbose is TRUE
#' @param ... parameters passed to cat
catv <- function(...) {
    if (verbose) {
        cat(...)
    }
}

#' Converts NA values in a vector to FALSE
#' @param x vector to convert NAs in
#' @return the x vector with NAs turned to FALSE
NA2FALSE <- function(x) {
    x[is.na(x)] <- FALSE;x
}

## Link to emptyDrops pre-print: 
## https://www.biorxiv.org/content/early/2018/04/04/234872

## Parse the input arguments
option_list <- list(
    make_option(c('-i','--input-rds'),
                type='character',
                default=NULL, ## required
                dest='input_rds',
                help='input RDS file containing the data matrix in dgCMatrix format gene x droplet orientation'),
    make_option(c('-o','--output-csv'),
                type='character',
                default=NULL, ## required
                dest='output_csv',
                help='output CSV file, if it exists the program will abort'),
    make_option(c('--fdr-cutoff'),
                type='numeric',
                default=0.10, ## FDR recommendation from paper
                dest='fdr_cutoff',
                help='FDR value cutoff for calling cells'),
    make_option(c('-v','--verbose'),
                type='logical',
                action='store_true',
                dest='verbose', 
                help='print verbose messages',
                default=FALSE),
    make_option(c('--transpose'),
		type='logical',
		dest='transpose',
		action='store_true',
		help='transpose the input matrix before processing',
		default=FALSE),
    make_option(c('--emptydrops-lower'),
                default=100, ## as per emptyDrops package
                help='emptydrops lower parameter',
                dest='ed_lower'),
    make_option(c('--emptydrops-niters'),
                default=10000, ## as per emptyDrops package
                help='emptydrops niters paramters',
                dest='ed_niters'),
    make_option(c('--emptydrops-test-ambient'),
                dest='ed_test_ambient',
                default=FALSE, ## as per emptyDrops package
                action='store_true',
                help='Empty drops parameter test.ambient'),
    make_option(c('--emptydrops-ignore'),
                default=NULL, ## as per emptyDrops package
                help='emptyDrops ignore parameter',
                dest='ed_ignore'),
    make_option(c('--emptydrops-alpha'),
                default=NULL, ## as per emptyDrops package
                help='emptyDrops alpha parameter',
                dest='ed_alpha'),
    make_option(c('--min-molecules'),
		default=1000,
		help='minimum number of molecules for a droplet to be called a cell',
		dest='min_molecules')
	
    ## TODO: Expose parallel functionality
)

## Parse the arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Check the parsed arguments
if(is.null(opt$input_rds))  errorExit("Input RDS is not specified\n")
if(is.null(opt$output_csv)) errorExit("Output CSV is not specified\n")
if(!file.exists(opt$input_rds)) errorExit("Input RDS doesn't exist!\n")
if(file.exists(opt$output_csv)) errorExit("Output CSV file exists!\n")
if(is.null(opt$min_molecules)) errorExit("Minimum number of molecules is not specified\n")

## Load the required libraries here
## NOTE: We do this after parsing arguments so that --help returns immediatedly
pkgs <- c('DropletUtils','Matrix', 'BiocParallel')
missing_pkgs <- pkgs[!pkgs %in% installed.packages()]
if(length(missing_pkgs) > 0) {
    cat('Error the following package(s) are missing:', paste0(missing_pkgs, collapse=', '),'\n')
    quit(save="no",status=1)
}

## Load the libs, suppressing warnings and other messages
suppressWarnings(suppressPackageStartupMessages(library('DropletUtils',quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library('Matrix',quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library('BiocParallel',quietly=TRUE)))

## Set script parameters, this need to be done after package loading
## as SerialParam() is not otherwise availbale
inputRDS <- opt$input_rds
output_csv <- opt$output_csv
FDRcutoff <- opt$fdr_cutoff
verbose <- opt$verbose
ed_param_lower <- opt$ed_lower
ed_param_niters <- opt$ed_niters
ed_param_test.ambient <- opt$ed_test_ambient
ed_param_ignore <- opt$ed_ignore
ed_param_alpha <- opt$ed_alpha
ed_param_BPPARAM <- SerialParam() ## TODO: Allow serial or parallele, with ncore specification
min_molecules <- opt$min_molecules

## Read the input file
inputMatrix <- readRDS(inputRDS)

## Check the input matrix
## Check class
if(!class(inputMatrix) %in% c( 'dgCMatrix','dgRMatrix' )) {
    cat(paste0('Error: input matrix is not of class dgCMatrix or dgRMatrix. It is of class:',class(inputMatrix)),file=stderr())
    quit(save="no",status=1)
} else {
    catv('Note: Input matrix class check: OK')
}

## Check dimensions
if(any(dim(inputMatrix) == 0)) {
    ## If the matrix is empty we can't run emptyDrops, write a empty table header
    cat('Warning: one or more dimensions of the input matrix are empty. Generating empty result table.',file=stderr())
    outputFile <- file(output_csv)
    writeLines('"CellId","Total","LogProb","PValue","Limited","FDR","IsCell"',outputFile)
    close(outputFile)
    q(save="no",status=0)
} else {
    catv('Note: Input matrix input dimensions check: OK')
}

## Check row and column names
if(is.null(rownames(inputMatrix))) {
    cat('Warning: rownames of the input matrix are empty',file=stderr())
    quit(save="no",status=1)
} else {
    catv('Note: Input matrix rownames check: OK')
}
if(is.null(colnames(inputMatrix))) {
    cat('Warning: colnames of the input matrix are empty',file=stderr())
    quit(save="no",status=1)
} else {
    catv('Note: Input matrix colnames check: OK')	
}

## Convert to CsparseMatrix if required
## emptyDrops seems to work fine with dgRMatrix so this is not required for now
## however it is here as it can help optimization in the future
## if (class(inputMatrix) == 'dgRMatrix') {
##   catv('Input Matrix is in dgRMatrix format. Converting...')
##   inputMatrix <- as(inputMatrix, "CsparseMatrix")
##   catv('done\n')
## }

## If requested transpose the input matrix
## NOTE: dgRMatrix to dgCMatrix conversion and transposition
## can be done in a single step by re-interpreting the indexes
if (opt$transpose) {
   catv('Transposing input matrix...')
   inputMatrix <- Matrix::t(inputMatrix)
   catv('done\n')
}

## Run emptyDrops with error handling
catv('Running emptyDrops...')
t0 <- Sys.time()
tryCatch({
    emptyDrops_result <- emptyDrops(m=inputMatrix,
                                    lower=ed_param_lower,
                                    niters=ed_param_niters,
                                    test.ambient=ed_param_test.ambient,
                                    ignore=ed_param_ignore,
                                    alpha=ed_param_alpha,
                                    BPPARAM=ed_param_BPPARAM)
},error=function(e) {
    if(grepl("need at least four unique 'x' values",e$message,fixed=TRUE))
    {
        cat('Error: an error occured while running emptyDrops!\n',file=stderr())
        cat('Error: ', e$message,'\n',file=stderr())
        # Write an empty_drops_results.csv file that has only NAs
        n_rows = dim(inputMatrix)[2] # Get number of rows
        emptyDrops_result <- matrix(data=NA,nrow=n_rows,ncol=7)
        ## Convert output from DataFrame to data.frame
        emptyDrops_result <- as.data.frame(emptyDrops_result)
        colnames(emptyDrops_result) <-  c("CellId","Total", "LogProb", "PValue", "Limited", "FDR", "IsCell")
        emptyDrops_result[,"CellId"] = colnames(inputMatrix)
        ## Write the output file
        catv('Writing output CSV with NA\'s instead of emptydrops metrics ...')
        write.csv(x=emptyDrops_result, file=output_csv,row.names=FALSE)
        catv('done\n')
        quit(save="no",status=0)
    }
    else
    {
        cat('Error: an error occured while running emptyDrops!\n',file=stderr())
        cat('Error: ', e$message,'\n',file=stderr())
        quit(save="no",status=1)
    }
})
t1 <- Sys.time()
emptyDrop_runtime <- t1 - t0
catv('done in', as.numeric(t1 - t0,units='secs'), 'seconds\n')

catv('Preparing output table...')
## Convert output from DataFrame to data.frame
emptyDrops_result <- as.data.frame(emptyDrops_result)

## Add the cell identifier as a column and remove from rownames
emptyDrops_result$CellId <- rownames(emptyDrops_result)
rownames(emptyDrops_result) <- NULL

## Call the cells according to the cutoff
emptyDrops_result$IsCell <- NA2FALSE(emptyDrops_result$FDR < FDRcutoff & emptyDrops_result$Total >= min_molecules)

## Define a column order
colOrder <- c("CellId","Total", "LogProb", "PValue", "Limited", "FDR", "IsCell")
emptyDrops_result <- emptyDrops_result[,colOrder]
catv('done\n') # Preparing output matrix

## Write the output file
catv('Writing output CSV...')
write.csv(x=emptyDrops_result, file=output_csv,row.names=FALSE)
catv('done\n')
