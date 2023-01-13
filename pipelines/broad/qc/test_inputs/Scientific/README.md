# A description of (some of) the test files here.
- ArrayVCFDownloadFP
    - Uses an Arrays VCF (hg19)
    - Downloads the fingerprint for the sample from Mercury.
- ArrayVCFSupplyFP
    - Uses an Arrays VCF (hg19)
    - Doesn't pull the FP from Mercury, it is supplied as an input to the test
- ArrayVCFSupplyFP
    - Uses an Arrays VCF (hg19)
    - Doesn't pull the FP from Mercury, it is supplied as an input to the test
    - Not supplying a fingerprint index file to test functionality
- JustDownloadFP_HG19
    - Only download the fingerprint genotypes from Mercury.
    - Output file on hg19
- JustDownloadFP_HG38
    - Only download the fingerprint genotypes from Mercury.
    - Output file on hg38
- ExomeCramDownloadFP_HG38
    - Uses a (big) CRAM from our Exome Scientific test (hg38)
    - Downloads the fingerprint for the sample from Mercury.
- ExomePlumbingCramDownloadFP_HG38
    - Uses a (small/downsampled) CRAM from our Exome plumbing test (hg38)
    - Doesn't pull the FP from Mercury, supply an NA2878 fingerprint
    - NOTE: Since this is a downsampled CRAM, the fingerprint LOD is poor (-1.8)
- PassThroughFP
    - Doesn't pull the FP from Mercury, just passes the input fingerprint VCF to output.
    - Just because.
- TcapBamDownloadFP_HG38
    - Uses a BAM from the TCAP pipeline
    - Downloads the fingerprint for the sample from Mercury.
- TcapBamSupplyFP_HG38
    - Uses a BAM from the TCAP pipeline
    - Doesn't pull the FP from Mercury, supply a fingerprint
- TcapPlumbingBamDownloadFP_HG38
    - Uses a (small/downsampled) BAM from the TCAP pipeline
    - Downloads the fingerprint for the sample from Mercury.
- TcapPlumbingBamSupplyFP_HG38
    - Uses a (small/downsampled) BAM from the TCAP pipeline
    - Doesn't pull the FP from Mercury, supply a fingerprint

    


