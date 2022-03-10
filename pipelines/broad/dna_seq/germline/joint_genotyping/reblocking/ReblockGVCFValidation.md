Some bugs in previous reblocking "truth", including phasing swap, 
so compare with germline single sample output.

Run ValidateVariants from GATK branch ldg_checkGTs:
- Reblocked output should have no gaps and no overlaps, i.e. should cover every 
position like a good GVCF 
- Reblocked output should not have two entries for any position
- Must contain all required annotations for GnarlyGenotyper and VQSR
-- some low coverage hets may be missing ReadPosRankSum or MQRankSum

Reblocked output will differ from input in the following ways
- Alleles not included in the called genotype are dropped
- MIN_DP annotation is dropped for reference blocks
- Hom-ref blocks are now combined into bands [0, 10), [10, 20), [20, 30), [30, 40), [40, 99]
- Hom-ref GQs are now "floored" such that the GQ of each block is the minimum for the appropriate band
- PLs are omitted for hom-ref genotypes
- Each position should be covered exactly once (except overlapping variants)

Variants GQs and likelihoods will stay the same

Extensive comparison code in the GATK branch ldg_VCFComparator (c4a1447ff9796e2bd48c2e3aafca3a50969a9749):
- Reblocking may drop variants, but won't introduce new variant positions
- HaplotypeCaller (HC) output variants that are not high quality (contain a * or <NON_REF> in genotype,
 are homozygous reference, have QUAL~0, are missing PLs) do not need to have a match in the reblocked output
- Genotype alleles should match
- HC output variants that are overlapped by a hom ref deletion can be dropped
- If there is no overlapping deletion for the HC variant, reblocked genotype does not need a *
- HC genotype alleles can have extra base padding, but trimmed representation should match
- Single base ref blocks now have END tags
 
  