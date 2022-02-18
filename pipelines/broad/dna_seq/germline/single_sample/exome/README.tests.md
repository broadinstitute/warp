### Germline Single Sample Workflow
## Exome Test Data

Output GVCFs will contain entries (incl. reference blocks) for all intervals in the 
exome calling intervals. For sample NA12878, reads likely derive from chr20 only, so
 expect a single GQ0 reference block per target over most of the exome with the first variant at chr7 (DP 3-5).
Good coverage spans roughly the interval chr20:10000000-31406318 with many reference blocks
and variant calls. A couple stray variants occur on chr22. 