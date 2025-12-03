# aou_9.0.1
2025-12-02

* This is the version of the pipeline used to create the ChrX VCF used for phasing in AoU v9 changes include:
    * Sanitize the chr name for the cluster creation. Cluster names may not include capital letters
    * Changed the cluster configuration from 2 dedicated workers and 50 preemptible workers to 25 dedicated and 25 preemptible workers
    * Increased time to live to 14400 min or 10 days. Actual run time was 69 hr 26 min or 4166 min. 

# aou_9.0.0
2025-10-02 (Date of Last Commit)

* This the version of the pipeline used to create the autosomal VCF files used for phasing in AoU v9 