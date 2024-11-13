## Announcement: CEMBA is Deprecated 9/12/2024

The CEMBA workflow is deprecated and no longer supported. File paths in the JSON file are also no longer supported.

However, the CEMBA documentation is still available. See [CEMBA Pipeline Overview](https://broadinstitute.github.io/warp/docs/Pipelines/CEMBA_MethylC_Seq_Pipeline/README) on the [WARP documentation site](https://broadinstitute.github.io/warp/)! The CEMBA data is also available on the NEMO data portal. The whitelist includes the following barcodes: CTCACG, CAGATC, CGATGT, ACTTGA, TTAGGC, GATCAG, TGACCA, TAGCTT, ACAGTG, GGCTAC, GCCAAT, CTTGTA.

### CEMBA summary

CEMBA is a pipeline developed by the [BRAIN Initiative](https://braininitiative.nih.gov/) that supports processing of multiplexed single-nuclei bisulfite sequencing data. It is an alignment and methylated base calling pipeline that trims adaptors, attaches cell barcodes, aligns reads to the genome, filters reads based on quality and creates a VCF with methylation-site coverage. 