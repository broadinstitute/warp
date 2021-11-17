---
sidebar_position: 14
---

# Reblocked GVCF Technical Overview
The following is a technical overview of the current reblocked GVCF schema (gnomADv4) as well descriptions of the previous GVCF schemas.

## gnomADv4 schema (4 bands)
* New “current” as of August 2021
* Arguments -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 

#### Summary (improvements compared with HC GVCFs):
* NOTE: this format is currently incompatible with GenotypeGVCFs!  Only GnarlyGenotyper can be specified in Joint genotyping WDL input json.
* PLs are omitted for homozygous reference sites to save space (which is why this format is currently incompatible with GenotypeGVCFs) -- GQs are output for genotypes, PLs can be approximated as [0, GQ, 2*GQ].

* Compared with old HaplotypeCaller GVCFs, GQ resolution for homozygous reference genotypes is reduced -- i.e. homRef GQs will be underconfident, which may affect analyses like de novo calling where confident reference genotypes are important.

* Alleles that aren’t called in the sample genotype are dropped, i.e. each variant should have no more than two alt alleles (in very rare cases, usually just one plus <NON_REF>)
New annotations to enable merging data for filtering without using genotypes -- RAW_GT_COUNT(S?) for doing ExcessHet calculation from a sites-only file; QUALapprox and/or AS_QUALapprox for doing QUAL approximation/filling QUAL VCF field from a combined sites-only field; VarDP and/or AS_VarDP used to calculate QualByDepth/QD annotation for VQSR
* No more MIN_DP -- joint DP varies?
* InbreedingCoeff values will vary due to some low quality genotypes being rounded down to GQ0

#### Cost/scale improvements

* Storage footprint is reduced compared with HaplotypeCaller output
* Fewer VariantContexts, i.e. lines, per VCF speeds up GenomicsDB/Hail import
* Fewer alternate alleles reduces memory requirements for merging


### Specific improvements compared with 7-band schema
* Not dropping GQ0s -- reblocked GVCF should cover all positions (that input GVCF covers)
* No overlaps -- only overlapping positions should be two variants, i.e. deletions, on separate haplotypes
(No more no-calls in GVCFs-- all genotypes should be called, positions with no data will be hom ref with GQ0; NOTE that this is not being validated and GenotypeGVCFs outputs will have no-calls to preserve previous behavior)
* Latest Terra WDL and config at https://app.terra.bio/#workspaces/ccdg-gatk-supp-terra/CCDG_WES_Reblocking_Work_F3/workflows/ccdg-gatk-supp-terra/ReblockGVCF-gatk4_exomes_goodCompression_fixed (NOTE that this has an off-by-one issue for contigs starts of exactly 1, which we don’t see in Broad production)

## CCDG WGS (7 bands)
* Suitable for de novo calling from exomes
* Arguments -drop-low-quals -do-qual-approx --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60
Q* UALapprox output in INFO field, required for GnarlyGenotyper
* Low quality alt alleles will be dropped (below stand-call-conf = 30)
* Alleles that aren’t called will be dropped 
* homRef calls at GQ0 (default) are dropped
* Non-ref AD zeroed out
* GTs with non-ref are converted to GQ0 homRef and dropped
* Blocks are floored -- GQ output for block is the lower bound for that block, regardless of median or min
* Current WDL at https://portal.firecloud.org/?return=terra#methods/methodsDev/ReblockGVCF-gatk4_exomes_goodCompression/4
* As used for UKBB


### gnomADv3 schema (deprecated -- 2 bands)
Arguments -drop-low-quals -do-qual-approx -rgq-threshold 10

* Default GQ blocks, < 20 and >= 20
* homRef calls below GQ10 are dropped
* → remaining tranche of homRef calls between GQ 10 and 20 gnomAD team had to filter
* QUALapprox output in INFO field, required for GnarlyGenotyper
* Low quality alt alleles will be dropped (below stand-call-conf = 30)
* Extra INFO annotations dropped, only data required for VQSR retained: DP, MQ, ReadPosRankSum, MQRankSum (strand bias data is in FORMAT)
* Non-ref AD zeroed out
* WDL at https://portal.firecloud.org/?return=terra#methods/eric-methods/ReblockGVCF-gatk4/1/wdl

## Original HaplotypeCaller schema (current Warp single-sample pipeline)
for (int i=1; i<=60; ++i) {
   GVCFGQBands.add(i);
}
GVCFGQBands.add(70); GVCFGQBands.add(80); GVCFGQBands.add(90); GVCFGQBands.add(99);
[1,2,3,4,5,....58,59,60, 70, 80, 90, 99]

