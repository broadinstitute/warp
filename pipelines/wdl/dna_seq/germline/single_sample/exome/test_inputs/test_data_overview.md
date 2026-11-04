# Exome Germline Single Sample Test Data
The following describes the data used for workflow validation, including workflow engineering (plumbing) tests and scientific validation tests.

* Output GVCFs will contain entries (incl. reference blocks) for all intervals in the exome calling intervals. 
* For sample NA12878, reads likely derive from chr20 only, so expect a single GQ0 reference block per target over most of the exome with the first variant at chr7 (DP 3-5).
* Good coverage spans roughly the interval chr20:10000000-31406318 with many reference blocks
and variant calls. A couple stray variants occur on chr22. 

## Plumbing
* RP-929.NA12878

## Scientific
* RP-929.NA12878
    * RP-1535.NA17-308 is a copy of NA12878 that has ~25% contamination
* SynDip
* NA19238 Yoruba daughter
* CEU parents (NA12891 and NA12892)