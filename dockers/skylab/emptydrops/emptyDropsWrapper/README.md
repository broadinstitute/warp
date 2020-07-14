# emptyDropsWrapper

Reads in an rds file containing the count matrix in dgCMatrix format in gene x droplet orientation and produces a CSV file indicating cell filtering and other cell metadat as returned by the emptyDrops function.

## Requirements

This script requires the following R packages to be installed: 

* optparse
* DropletUtils

DropletUtils are available from bioconductor and require R version 3.5

https://bioconductor.org/packages/release/bioc/html/DropletUtils.html

## Testing

To run the tests run the following:

```
cd /tools/emptyDropsWrapper/test/ 
./test_emptyDropsWrapper.sh
```
