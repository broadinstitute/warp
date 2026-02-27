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
| UniqueAlignFinalReads | Final unique mapped total reads after picard deduplication. | 
| UniqueAlignDuplicatedReads | Paired and unpaired duplicated reads. |
| UniqueAlignPCRDuplicationRate | FinalReads /(FinalReads and DuplicatedReads). |
| CisContacts | Number of chromatin contacts where the two loci are on the same chromosome. |
| CisCutContacts | Number of read pairs that are split from the same read at the cut site, and map to the same chromosome. |
| CisMultiContacts | CisContacts read pair contains multiple read contacts. |
| CisCutMultiContacts | CisCutContacts read pair contains multiple read contacts. |
| TransContacts | Number of chromatin contacts where the two loci are on different chromosome. |
| TransCutContacts | Number of read pairs that are split from the same read at the cut site, and map to different chromosomes. |
| TransMultiContacts | TransContacts read pair contains multiple read contacts. |
| TransCutMultiContacts | TransCutContacts read pair contains multiple read contacts. | 
| ChimericContacts | Two reads that are split from the same read, but not at the cut site, this might be due to artificial chimeric synthesis event. |
| NoContacts | Not a contact. |
| MappedFragments | Total number of mapped fragments. |
| DeduppedContacts | Total number of deduplicated contacts. |
| ContactsDeduplicationRate | `(Input_contacts - dedup_contacts) / (input_contacts + 0.00001)` |
| TotalCisContacts | Total number of cis contacts. |
| TotalTransContacts | Total number of trans contacts. |
| TotalMultiContacts | Total number of multi contacts (read pair contains multiple read contacts). |
| CisContactsRatio | TotalCisContacts / number of mapped fragments. |
| TransContactsRatio | TotalTransContacts / number of mapped fragments. |
| MultiContactsRatio | TotalMultiContacts / No. of mapped fragments. |
| mCCCmC | Total methylated cytosine in the CCC context. |
| mCGmC | Total methylated cytosine in the CG context. |
| mCHmC | Total methylated cytosine in the CH context. |
| mCCCCov | Total covered cytosine in the CCC context. |
| mCGCov | Total covered cytosine in the CG context. |
| mCHCov | Total covered cytosine in the CH context. |
| mCCCFrac | Fraction of methylated cytosine (`mCCCmC`) divided by covered cytosine (`mCCCCov`) in the CCC context. |
| mCGFrac | Fraction of methylated cytosine (`mCGmC`) divided by covered cytosine (`mCGCov`) in the CG context. |
| mCHFrac | Fraction of methylated cytosine (`mCHmC`) divided by covered cytosine (`mCHCov`) in the CH context. |


