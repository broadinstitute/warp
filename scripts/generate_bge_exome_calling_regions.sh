#!/bin/bash

# Clinvar code below was adapted from Ben Weisburd in 2019
hg38=/seq/references/Homo_sapiens_assembly38/v0/
picard=/seq/software/picard/1.1510/bin/picard-private-05132dd410be897226a14eeba9fa47c8dd58d1cf-all.jar
gatk=/seq/software/gatk/gatk-4.0.5.1/gatk


# Get Twist targets from https://www.twistbioscience.com/sites/default/files/resources/2021-10/Twist%20Alliance%20Clinical%20Research%20Exome_Covered_Targets_hg38-34.9MB%20.bed
wget -nc https://www.twistbioscience.com/sites/default/files/resources/2021-10/Twist%20Alliance%20Clinical%20Research%20Exome_Covered_Targets_hg38-34.9MB%20.bed \
  -O TwistAllianceClinicalResearchExome_Covered_Targets_hg38.bed
java -jar $picard BedToIntervalList -I TwistAllianceClinicalResearchExome_Covered_Targets_hg38.bed \
  -SD $hg38/Homo_sapiens_assembly38.dict -O TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list
#Wrote 201834 intervals spanning a total of 34883444 bases

# ====================================================
# =============== gencode.v28.interval_list ==========
# ====================================================
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

# This should include CDS, utr, miRNA, lincRNA and wierd "antisense" exons
cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list
zcat gencode.v28.annotation.gtf.gz | awk '$3 == "exon"  { print $1, $4, $5, $7, $18}' | tr ' ' '\t' | tr -d '";' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    O=gencode.v28.interval_list \
    UNIQUE=true
rm tmp.interval_list
#Produced 310940 intervals totalling 133277685 bases.

# ==========================================================
# =============== mirbase.v22.interval_list ================
# ==========================================================

wget -nc -O mirbase.v22.gff https://www.mirbase.org/ftp/22/genomes/hsa.gff3

cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list
cat mirbase.v22.gff | awk '$3 == "miRNA_primary_transcript" {sub(/.*Name=/, "", $9); print $1,$4,$5,$7,$9}' | tr ' ' '\t' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    O=mirbase.v22.interval_list \
    UNIQUE=true
rm tmp.interval_list
#Produced 1859 intervals totalling 152866 bases.

# ==========================================================
# ============= illumina_exome.interval_list ===============
# ==========================================================
ice_targets_lifted=/seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list
cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list

grep 'blat' $ice_targets_lifted > blats.intervals
python removeBadBlats.py blats.intervals good.intervals
cp $hg38/Homo_sapiens_assembly38.dict header
cat header good.intervals > goodBlats.interval_list
egrep -v '^@' $ice_targets_lifted | egrep -v 'blat' \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, "ice_target_" FNR}' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    I=goodBlats.interval_list \
    O=ice_coding_v1_targets.interval_list \
    ACTION=UNION
rm tmp.interval_list
#Produced 216560 intervals totalling 37864250 bases.


wget -nc ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230121.vcf.gz
#remove variants that are not on the canonical autosomes and filter for just pathogenic
zgrep -i '^#\|pathogenic' clinvar_20230121.vcf.gz | grep -v '^MT' | sed $'s/^/chr/' | sed $'s/chr[#]/#/' > clinvar_20230121.GRCh38.pathogenic.vcf

$gatk UpdateVcfSequenceDictionary -I clinvar_20230121.GRCh38.pathogenic.vcf -SD $hg38/Homo_sapiens_assembly38.dict -O clinvar_20230121.GRCh38.pathogenic.fixed.vcf
$gatk IndexFeatureFile -F clinvar_20230121.GRCh38.pathogenic.fixed.vcf

#Limit clinvar variants to only non-coding (outside of what is already covered by other intervals)
#Remove deletions longer than 500 bases since we won't be able to call them from BGE short reads
$gatk SelectVariants -XL TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list \
  -XL ice_coding_v1_targets.interval_list -XL mirbase.v22.interval_list -XL gencode.v28.interval_list \
	-V clinvar_20230121.GRCh38.pathogenic.fixed.vcf -O clinvar_20230121_noncoding_no_long_dels_pathogenic.vcf \
	-select '!(vc.isIndel()&&vc.isBiallelic()&&vc.getIndelLengths().get(0)<-500)'

java -jar $picard IntervalListTools -I clinvar_20230121_noncoding_no_long_dels_pathogenic.vcf \
	-O clinvar_20230121_noncoding_non_long_deletion_pathogenic.interval_list
#Produced 16299 intervals totalling 164364 bases

#Add names to twist and clinvar interval lists, others already has this
#Name the twist intervals with an incremented number so they are unique
grep -v '^@' TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list | awk 'sub(/\./,"twist_target_"c++)' > TwistNames.intervals
grep '^@' TwistAllianceClinicalResearchExome_Covered_Targets_hg38.interval_list > twistHeader.txt
cat twistHeader.txt TwistNames.intervals > TwistAllianceClinicalResearchExome_Covered_Targets_hg38_named.interval_list

#Clinvar has a number in the name column (5th column) so leave it as is but add a ClinVar string
grep -v '^@' clinvar_20230121_noncoding_non_long_deletion_pathogenic.interval_list | awk 'BEGIN{OFS="\t"}$5="clinvar_pathogenic_"$5' > clinvarNames.intervals
grep '^@' clinvar_20230121_noncoding_non_long_deletion_pathogenic.interval_list > clinvarHeader.txt
cat clinvarHeader.txt clinvarNames.intervals > clinvar_20230121_noncoding_non_long_deletion_pathogenic_named.interval_list



#merge all five lists
java -jar $picard IntervalListTools -UNIQUE \
	-I gencode.v28.interval_list \
	-I mirbase.v22.interval_list \
	-I ice_coding_v1_targets.interval_list \
	-I clinvar_20230121_noncoding_non_long_deletion_pathogenic_named.interval_list \
	-I TwistAllianceClinicalResearchExome_Covered_Targets_hg38_named.interval_list \
	-PADDING 50 \
	-O merged_all_lists.interval_list

#Produced 296381 intervals totalling 165239993 bases.

#Limit intervals to chr1-22,chrX,chrY (remove alts and chrM)
echo "chrY	1	57227415	+	." > chrY.tmp
cat /seq/references/Homo_sapiens_assembly38/v0/resources/wholegenome.interval_list chrY.tmp > main_wholegenome.interval_list

java -jar $picard IntervalListTools -ACTION INTERSECT \
	-I merged_all_lists.interval_list \
	-I main_wholegenome.interval_list \
	-O bge_exome_calling_regions.v1.1.interval_list

#Produced 296328 intervals totalling 165207924 bases

#convert to bed for Dragen
java -jar $picard IntervalListToBed -I bge_exome_calling_regions.v1.1.interval_list \
  -O bge_exome_calling_regions.v1.1.bed
