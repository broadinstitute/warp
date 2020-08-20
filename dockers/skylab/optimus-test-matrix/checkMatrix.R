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


## Here we are checking the matrices for equality by looking at the
## element contents which will be identical after drop0() has been run for
## identical matrices
if(all(matrix@i == referenceMatrix@i) && all(matrix@p == referenceMatrix@p) &&
  all(matrix@Dimnames[[1]] == referenceMatrix@Dimnames[[1]]) &&
  all(matrix@Dimnames[[2]] == referenceMatrix@Dimnames[[2]]) &&
  all(matrix@x == referenceMatrix@x)) {
  print('PASS: Matrices are identical')
  quit(status=0)
} else {
  print('FAIL: Matrices differ')
  quit(status=1)
}