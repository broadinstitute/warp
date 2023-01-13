#!/usr/bin/env Rscript

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Date: December 12, 2018
## Description: An R script that will read text converted npz/npy arrays
##   from the optimus pipeline and save it as an R rds file that contains a
##   single dgRMatrix matrix (from the Matrix package)

## Parse the input arguments
library('optparse')
option_list = list(
    make_option(c('-i','--input-dir'), type='character', default=NULL,help='input directory'),
    make_option(c('-o','--output-file'), type='character',default=NULL,help='output rds file')
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check arguments
if (is.null(opt$`input-dir`)) stop('input directory is not specified');
if (is.null(opt$`output-file`)) stop('output file is not specified');

inputdir <- opt$`input-dir`
outfile <- opt$`output-file`

#################################
## Main Script
library(Matrix)

col_index_file <- file.path(inputdir, 'sparse_counts_col_index.txt')
indices_file <- file.path(inputdir, 'sparse_counts_indices.txt')
row_index_file <- file.path(inputdir, 'sparse_counts_row_index.txt')
data_file <- file.path(inputdir, 'sparse_counts_data.txt')
indptr_file <- file.path(inputdir, 'sparse_counts_indptr.txt')
shape_file <- file.path(inputdir, 'sparse_counts_shape.txt')



#' Get the filesize of the specified 
getFileSize <- function(filename) {
    file.info(filename)$size
}

## Flag for empty matrix -- no need to sort if empty
matrixIsEmpty <- FALSE

cat('Loading text files...')
tryCatch({
    shape_file_size <- getFileSize(shape_file)
    if(shape_file_size >0) {
        shape <- read.table(shape_file, stringsAsFactors=FALSE)$V1
    } else {
        cat('\nAn error occured while reading the input shape file. Exiting.');
        ## This is fatal, exit with error code
        q(status=1);
    }

    if (all(shape == 0)) {
        ## This is an empty matrix, warn and continue
        cat('\nWarning: Loading an empty matrix!\n')
        matrixIsEmpty <- TRUE
        col_index <- c()
        row_index <- c()
        data <- double()
        indices <- integer()
        indptr <- as.integer(c(0))
    } else {
        col_index <- read.table(col_index_file, stringsAsFactors=FALSE)$V1
        row_index <- read.table(row_index_file, stringsAsFactors=FALSE)$V1
        data <- read.table(data_file, stringsAsFactors=FALSE)$V1
        indices <- read.table(indices_file, stringsAsFactors=FALSE)$V1
        indptr <- read.table(indptr_file, stringsAsFactors=FALSE)$V1
    }
    
}, error=function(e) {
    cat(paste0("\nAn error occured while loading the files: ",e," Exiting.\n"));
    q(status=1);
}) 
cat('done\n')

## The R Matrix package requires that the second dimension is ordered
## within the first dimention, but the input isn't, here we output a sort order
## by sorting all the second index dimentions for each range of first index
if(!matrixIsEmpty) {
    cat('Sorting elements...')
    ord1 <- unlist(lapply(1:(length(indptr)-1), function(i) {
        rowElementOrder<-vector(mode="list", length=0);
        if (indptr[i]!=indptr[i+1]) {
           rowStart <- indptr[i] + 1;
           rowEnd <- indptr[i+1];
           rowIndices <- indices[c(rowStart:rowEnd)]
           rowElementOrder <- order(rowIndices) + rowStart - 1;
        }
        rowElementOrder
    }))
    cat('done\n')
} else {
    ord1 <- c()
}

dimError=FALSE;
if ( shape[1] != length(row_index) ) {
    cat('Error: the number of rows in the matrix does not equal the number of rows in the index.\n');
    cat("    Number of array entries is ",shape[1],", number of labels",length(row_index), "\n");
    dimError=TRUE;
}

if ( shape[2] != length(col_index) ) {
    cat('Error: the number of rows in the matrix does not equal the number of rows in the index.\n');
    cat("    Number of array entries is ",shape[2],", number of labels",length(col_index), "\n");
    dimError=TRUE;
}

if (dimError) {
    stop('Dimension error');
}

cat('Lengths: ', 'indices=', length(indices),  '  ord1 = ', length(ord1), ' data = ', length(data),
        ' indptr = ', length(indptr), 
        ' shape = ', shape, ' row index = ', length(row_index), ' col index = ', length(col_index), '\n')
## Reorder the matrix (j and x values as per the ordering done above)
## and generate the dgRMatrix object
cat('Generating matrix...')
res <- new('dgRMatrix',
           j=as.integer(indices[ord1]),
           p=as.integer(indptr),
           x=as.double(data[ord1]),
           Dim=as.integer(shape),
           Dimnames=list(row_index,col_index)
           )
cat('done\n')

cat('Saving matrix...')
saveRDS(object=res,file=outfile)
cat('done\n')

### End of Script
