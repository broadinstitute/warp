#!/usr/bin/env R

## Load Matrix Library
library(Matrix)

## Read the matrix
cat('Reading input matrix...')
matrix <- readRDS('matrix.rds')
cat('done\n')

## Print info on matrix dimensions
matrix.dims <- dim(matrix)
print(paste0('Input matrix dimensions: ', matrix.dims[1] , ' x ', matrix.dims[2],'\n'))

## Put rows and columns in a defined order
cat('Ordering matrix columns and rows...')
matrix <- matrix[ sort(rownames(matrix)), sort(colnames(matrix)) ]
matrix <- drop0(matrix)
cat('done\n')

## We generate a new reference matrix here in case we want to replace
cat('Generating new reference matrix...')
saveRDS(matrix,'newReferenceMatrix.rds')
gc()
cat('done\n')

## Generate simple diagnostic plots
cat('Generating diagnostic plots...')
png('reads_per_cell_histogram.png')
hist(rowSums(matrix),main="Reads per Cell")
dev.off()

png('reads_per_gene_histogram.png')
hist(colSums(matrix), main="Reads per Gene")
dev.off()

png('number_of_genes_per_cell.png')
hist(colSums(matrix > 1), main="Genes per Cell")
dev.off()
cat('done\n')

## Read in the reference matrix
cat('Reading in reference matrix...')
referenceMatrix <- readRDS('referenceMatrix.rds')
cat('done\n')

cat('Ordering reference matrix columns and rows...')
referenceMatrix <- referenceMatrix[ sort(rownames(referenceMatrix)), sort(colnames(referenceMatrix)) ]
cat('done\n')

ref.matrix.dims <- dim(referenceMatrix)
cat(paste0('Input matrix dimensions: ', ref.matrix.dims[1] , ' x ', ref.matrix.dims[2],'\n'))

## Check if the matrices are identical by looking at the element contents
## Give specific error messages for all the differences before quitting
different <- FALSE

if(!all(matrix@i == referenceMatrix@i)){
  cat('FAIL: i differs between the two matrices\n')
  different <- TRUE
}

if(!all(matrix@p == referenceMatrix@p)){
  cat('FAIL: p differs between the two matrices\n')
  different <- TRUE
}

if(!all(matrix@Dimnames[[1]] == referenceMatrix@Dimnames[[1]])){
  cat('FAIL: Dimnames[[1]] differs between the two matrices\n')
  different <- TRUE
}

if(!all(matrix@Dimnames[[2]] == referenceMatrix@Dimnames[[2]])){
  cat('FAIL: Dimnames[[2]] differs between the two matrices\n')
  different <- TRUE
}

## x can have 0.1% differences
if(sum((matrix@x-referenceMatrix@x)!=0)/length(matrix@x) * 100 > 0.1) {
  cat('FAIL: There are too many differences in x between the two matrices\n')
  different <- TRUE
}

if(different) {
  quit(status=1)
} else {
  cat('PASS: Matrices are identical')
  quit(status=0)
}
