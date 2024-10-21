# STAR Aligner Metrics
The STAR aligner produces multiple text files containing library-level summary metrics, cell-level metrics, and UMI metrics. The Optimus workflow compresses these files into a single TAR. These outputs are directly from the aligner as different batches of the data are analyzed in parallel. 

The STAR aligner metrics are supplemental to the [library-level metrics CSV](./Library-metrics.md) that is also produced by Optimus. Several of the calculations produced in the library metrics are directly based on the STAR aligner metrics. 

The following sections describe these outputs. 

## Align Features Metrics 
The Align feature text file contains library-level metrics produced by the STARsolo alignment detailing the alignment of reads to genomic features during single-cell RNA-seq analysis. These metrics indicate how well reads map to specific genomic features or whether they failed to map due to various reasons. For example: 
**noUnmapped** represents the number of reads that were not aligned to any feature in the genome.
**noNoFeature** reflects reads that were aligned but did not map to any specific feature such as exons or genes.
**MultiFeature** counts reads that were aligned to multiple features.
**yesWLmatch** and **yesCellBarcodes** track how well reads match the barcode whitelist, an important step in identifying valid cell barcodes, which helps demultiplex the single-cell RNA-seq data​.

Each of the table metrics gives insights into different stages of read alignment, from barcode matching to gene feature mapping, allowing you to assess the quality and accuracy of the alignment step in the pipeline.


| Metrics name | Description | 
| --- | --- | 
| noUnmapped | Number of unmapped reads | 
| noNoFeature | Number of reads not mapped to a feature. | 
| MultiFeature | Number of reads aligned to multiple features. | 
| subMultiFeatureMultiGenomic | Number of reads mapping to multiple genomic loci and multiple features. | 
| noTooManyWLmatches | Number of reads not counted because their barcoded pair has too many matches to the whitelist. | 
| noMMtoWLwithoutExact | Number of reads not counted because their barcoded pair has mismatches to the whitelist and there's no more reads supporting that barcode. | 
| yesWLmatch | Number of reads whose barcoded pair has a match to the whitelist. |
| yessubWLmatchExact | Number of reads with cell barcode exactly matched to the whitelist (a subset of yesWLmatch). | 
| yessubWLmatch_UniqueFeature | Number of reads matched to the WL and unique feature (a subset of yesWLmatch). | 
| yesCellBarcodes | Number of reads associated with a valid cell barcode. | 
| yesUMIs | Number of reads associated with a valid UMI. | 






## Cell Read Metrics 

The **cell read metrics** text file provides cell barcode-level information about the reads; for instance:
**cbMatch** counts the number of reads that successfully matched the cell barcode.
**cbPerfect** gives the number of reads with a perfect match to a cell barcode, while **cbMMunique** and **cbMMmultiple** measure mismatches that still align uniquely or to multiple barcodes, respectively.
**genomeU** and **genomeM** count reads mapped to one or multiple loci in the genome, respectively.
**exonic** and **intronic** track reads mapping to annotated exons or introns, helping distinguish between different gene regions in the analysis.

These metrics are important for assessing the quality of individual cell barcodes.

| Metrics  |  Description  | 
| --- | --- |        
| CB | Cell barcode  | 
| cbMatch | Number of reads that matched the cell barcode. |  
| cbPerfect | Number of perfect matches on cell barcode. |    
| cbMMunique | Number of reads with cell barcodes that map with mismatches to one barcode in the passlist. | 
| cbMMmultiple | Number of reads with cell barcodes that map with mismatches to multiple barcodes in the passlist. |    
| genomeU | Number of reads mapping to one locus in the genome. |  
| genomeM | Number of reads mapping to multiple loci in the genome. |  
| featureU | Number of reads mapping to one feature (Gene, GeneFull, etc). | 
| featureM | Number of reads mapping to multiple features. |  
| exonic | Number of reads mapping to annotated exons. | 
| intronic | Number of reads mapping to annotated introns; these are only calculated for --soloFeatures GeneFull_Ex50pAS and/or GeneFull_ExonOverIntron. | 
| exonicAS | Number of reads mapping antisense to annotated exons. | 
| intronicAS | Number of reads mapping antisense to annotated introns; these are only calculated for --soloFeatures GeneFull_Ex50pAS. | 
| mito | Number of reads mapping to the mitochondrial genome. |  
| countedU | Number of unique-gene reads whose UMIs contributed to counts in the matrix.mtx (eads with valid CB/UMI/gene). |        
| countedM | Number of multi-gene reads whose UMIs contributed to counts in the matrix.mtx. |    
| nUMIunique | Total number of counted UMI for unique-gene reads. |  
| nGenesUnique | Number of genes for unique-gene reads. |      
| nUMImulti | Total number of counted UMI for multi-gene reads. |  
| nGenesMulti | Number of genes for multi-gene reads. |        

## Summary.txt

The **summary** text file contains additional library-level metrics produced by the STARsolo aligner, such as:  
**Number of reads**, which reflects the total reads processed, and **reads with valid barcodes**, which indicates how many reads matched the barcode whitelist.
**Sequencing saturation** shows the completeness of sequencing, where higher values indicate fewer additional reads are needed to capture new UMIs.
Metrics like **Q30 Bases in CB+UMI** and **Q30 Bases in RNA read** give insights into sequencing quality, showing how many reads had high-quality base calls.
Other key metrics, such as **reads mapped to the genome: Unique+Multiple** and **estimated number of cells**, provide a sense of how well reads were mapped to the genome and how many cells were identified.
These summary metrics help users assess the overall quality and completeness of their single-cell RNA-seq data, serving as a useful checkpoint for determining whether the data is suitable for further analysis. 

| Metric | Description |
| --- | --- |
| Number of Reads | Number of reads in the library. |
| Reads With Valid Barcodes | Fraction of reads with valid barcodes. |
| Sequencing Saturation | Proportion of unique molecular identifiers (UMIs) that have been sequenced at least once compared to the total number of possible UMIs in the sample; calculated as: 1-(yesUMIs/yessubWLmatch_UniqueFeature). |
| Q30 Bases in CB+UMI | Fraction of high-quality reads in the cell barcode and UMI read.  |
| Q30 Bases in RNA read | Fraction of high-quality reads in the RNA read.  |
| Reads Mapped to Genome: Unique+Multiple | Fraction of unique and multimapped reads that mapped to the genome.  |
| Reads Mapped to Genome: Unique | Fraction of unique reads that mapped to the genome. |
| Reads Mapped to genes: Unique+Multiple | Fraction of reads that mapped to genes as defined by the –solo-feature parameter. |
| Reads Mapped to Genes: Unique| Fraction of unique reads that mapped to genes. |
| Estimated Number of Cells | Number of barcodes that STARsolo flagged as cells based on UMIs.  |
| Unique Reads in Cells Mapped to genes | Total number of unique reads that mapped to genes across all cells |
| Fraction of Unique Reads in Cells | Fraction of unique reads across all cells. |
| Mean Reads per Cell | Mean number of reads per cell. |
| Median Reads per Cell | Median number of reads per cell. |
| UMIs in Cells | Number of UMIs per cell. |
| Mean UMI per Cell | Mean number of UMIs per cell. |
| Median UMI per Cell | Median number of UMI per cell. |
| Mean Genes per Cell | Mean number of genes expressed per cell. |
| Median Genes per Cell | Median number of genes per cell. |
| Total Genes Detected | Total number of genes detected in the overall library. |


## UMI per cell
The UMI per cell text file is a list of UMI counts per every cell. It contains two columns. The first column contains the number of UMIs per each barcode entry. The second column indicates whether a barcode was flagged as a cell. A 1 indicates that it passed filtering criteria to be considered a cell and 0 indicates that it did not pass.
