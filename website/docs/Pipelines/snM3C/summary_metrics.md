---
sidebar_position: 2
---

# snm3C Mapping Summary Metrics Overview

The snm3C pipeline outputs a summary CSV file containing trimming, mapping, deduplication, chromatin contact, and ALLC site statistics. 

The summary file is generated using the `smn3c_summary()` function of a [custom python3 script](https://github.com/lhqing/cemba_data/blob/788e83cd66f3b556bdfacf3485bed9500d381f23/cemba_data/hisat3n/summary.py).

The snm3C pipeline was adapted from YAP (Yet Another Pipeline) in collaboration with Hanqing Liu, Wei Tian, Wubin Ding, Huaming Chen, Chongyuan Luo, Jingtian Zhou, and the entire laboratory of Joseph Ecker. For more information about the snm3C metrics, please see the [YAP documentation](https://hq-1.gitbook.io/mc/) created by Hanqing Liu.

| Metric | Details |
|:------ | :------ |
| cell | The unique identifier for each cell. |
| InputReadPairs | Total number of read pairs. |
| InputReadPairsBP | Total number of base pairs. |
| TrimmedReadPairs | Total number of trimmed read pairs. |
| R1WithAdapters | Number of R1 reads trimmed to remove Illumina adapter. |
| R1QualTrimBP | Number of R1 base pairs trimmed due to low base quality. |
| R1TrimmedReadsBP | Number of R1 base pairs remaining after adapter and quality trimming. |
| R2WithAdapters | Number of R2 reads trimmed to remove Illumina adapter. |
| R2QualTrimBP | Number of R2 base pairs trimmed due to low base quality. |
| R2TrimmedReadsBP | Number of R2 base pairs remaining after adapter and quality trimming. |
| UniqueMappedReads | Number of uniquely mapped reads. |
| UniqueMappingRate | Rate of unique read mapping. |
| MultiMappedReads | Number of multimapped reads. | 
| MultiMappingRate | Rate of multimapping. |
| OverallMappingRate | Rate of mapping for all reads. |
| R1SplitReadsUniqueMappedReads | Number of uniquely mapped R1 reads. |
| R1SplitReadsUniqueMappingRate | Rate of unique R1 read mapping. |
| R1SplitReadsMultiMappedReads | Number of multimapped R1 reads. | 
| R1SplitReadsMultiMappingRate | Rate of R1 read multimapping. |
| R1SplitReadsOverallMappingRate | Rate of mapping for all R1 reads. |
| R2SplitReadsUniqueMappedReads | Number of uniquely mapped R2 reads. |
| R2SplitReadsUniqueMappingRate | Rate of unique R2 read mapping. |
| R2SplitReadsMultiMappedReads | Number of multimapped R2 reads. | 
| R2SplitReadsMultiMappingRate | Rate of R2 read multimapping. |
| R2SplitReadsOverallMappingRate | Rate of mapping for all R2 reads. |
| UniqueAlignFinalReads | | 
| UniqueAlignDuplicatedReads | |
| UniqueAlignPCRDuplicationRate |
| CisContacts | Number of chromatin contacts where the two loci are on the same chromosome. |
| CisCutContacts | 
| CisMultiContacts | 
| CisCutMultiContacts |
| TransContacts | Number of chromatin contacts where the two loci are on different chromosome. |
| TransCutContacts |
| TransMultiContacts |
| TransCutMultiContacts |
| ChimericContacts |
| NoContacts |
| MappedFragments |
| DeduppedContacts |
| ContactsDeduplicationRate |
| TotalCisContacts |
| TotalTransContacts |
| TotalMultiContacts |
| CisContactsRatio |
| TransContactsRatio |
| MultiContactsRatio |
| mCCCmC | Total methylated cytosine in the CCC context. |
| mCGmC | Total methylated cytosine in the CG context. |
| mCHmC | Total methylated cytosine in the CH context. |
| mCCCCov | Total covered cytosine in the CCC context. |
| mCGCov | Total covered cytosine in the CG context. |
| mCHCov | Total covered cytosine in the CH context. |
| mCCCFrac | Fraction of methylated cytosine (`mCCCmC`) divided by covered cytosine (`mCCCCov`) in the CCC context. |
| mCGFrac | Fraction of methylated cytosine (`mCGmC`) divided by covered cytosine (`mCGCov`) in the CG context. |
| mCHFrac | Fraction of methylated cytosine (`mCHmC`) divided by covered cytosine (`mCHCov`) in the CH context. |


