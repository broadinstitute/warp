Some bugs in previous reblocking "truth", including phasing swap, 
so compare with germline single sample output.

Run ValidateVariants from GATK branch ldg_checkGTs:
- Reblocked output should have no gaps and no overlaps, i.e. should cover every 
position like a good GVCF 
- Reblocked output should not have two entries for any position
- Must contain all required annotations for GnarlyGenotyper and VQSR
-- some low coverage hets may be missing ReadPosRankSum or MQRankSum

Reblocked output will differ from input

Extensive comparison code in the GATK branch ldg_VCFComparator:
- Reblocking may drop variants, but won't introduce new variant positions
- HaplotypeCaller (HC) output variants that are not high quality (contain a * or <NON_REF> in genotype,
 are homozygous reference, have QUAL~0, are missing PLs) do not need to have a match in the reblocked output
- Genotype alleles should match
- HC output variants that are overlapped by a hom ref deletion can be dropped
- If there is no overlapping deletion for the HC variant, reblocked genotype does not need a *
- HC genotype alleles can have extra base padding, but trimmed representation should match
- Single base ref blocks now have END tags
 
  