###############################################################################
# README file for how the calling regions for the single-sample HC based
# exome calling pipeline were created.  The following is a two-part process
# of identifying and formating various input sets of calling regions and then
# secondarily combining and padding said regions.
###############################################################################

hg38=/seq/references/Homo_sapiens_assembly38/v0/
picard=/seq/software/picard/1.1368/bin/picard-private-d5b94667a33ea4b110def476023b994de25868ef-all.jar

# ====================================================
# =============== gencode.v28.interval_list ==========
# ====================================================
if [ ! -e gencode.v28.annotation.gtf.gz ]; then 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
fi

# This should include CDS, utr, miRNA, lincRNA and wierd "antisense" exons
cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list
zcat gencode.v28.annotation.gtf.gz | awk '$3 == "exon"  { print $1, $4, $5, $7, $18}' | tr ' ' '\t' | tr -d '";' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    O=gencode.v28.interval_list \
    UNIQUE=true
rm tmp.interval_list


# ==========================================================
# =============== mirbase.v22.interval_list ================
# ==========================================================

if [ ! -e mirbase.v22.gff ]; then 
wget -O mirbase.v22.gff ftp://mirbase.org/pub/mirbase/22/genomes/hsa.gff3
fi

cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list
cat mirbase.v22.gff | awk '$3 == "miRNA_primary_transcript" {sub(/.*Name=/, "", $9); print $1,$4,$5,$7,$9}' | tr ' ' '\t' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    O=mirbase.v22.interval_list \
    UNIQUE=true
rm tmp.interval_list


# ==========================================================
# ============= illumina_exome.interval_list ===============
# ==========================================================
ice_targets_lifted=/seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list \
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

# ===============================================
# ======= Put it all together and pad it ========
# ===============================================
java -jar $picard IntervalListTools \
     I=gencode.v28.interval_list \
     I=mirbase.v22.interval_list \
     I=ice_coding_v1_targets.interval_list \
UNIQUE=true \
     PADDING=50 \
     O=exome_calling_regions.v1.1.interval_list


# ====================================================
# =============== gencode.v28.interval_list ==========
# ====================================================
if [ ! -e gencode.v28.annotation.gtf.gz ]; then 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
fi

cp $hg38/Homo_sapiens_assembly38.dict tmp.interval_list
zcat gencode.v28.annotation.gtf.gz | awk '$3 == "CDS" && ($26=="1;" || $26=="2;")  { print $1, $4, $5, $7, $18}' | tr ' ' '\t' | tr -d '";' >> tmp.interval_list
java -jar $picard IntervalListTools \
    I=tmp.interval_list \
    O=exome_evaluation_regions.v1.interval_list \
    UNIQUE=true
rm tmp.interval_list

