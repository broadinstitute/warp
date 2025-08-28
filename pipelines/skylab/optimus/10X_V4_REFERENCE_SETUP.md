# 10x 3' v4 (GEM-X) Chemistry Support for Optimus

This document explains how to set up and use the new 10x 3' v4 (GEM-X) chemistry support in Optimus.

## What's New

Optimus now supports 10x 3' v4 (GEM-X) chemistry in addition to the existing v2 and v3 chemistries. The v4 chemistry uses:
- 16-bp cell barcodes (same as v2/v3)
- 12-bp UMIs (same as v3)
- New barcode inclusion list: `3M-3pgex-may-2023.txt.gz`

## Required Downloads

### 1. Barcode Inclusion List (Whitelist)

**File**: `3M-3pgex-may-2023.txt.gz` or `3M-3pgex-may-2023_TRU.txt.gz`

Aug 28, 2025: These files are temporarily available here:
 
- gs://new-files-for-reference/3M-3pgex-may-2023.txt.gz
- gs://new-files-for-reference/3M-3pgex-may-2023_TRU.txt.gz
- gs://new-files-for-reference/3M-5pgex-jan-2023.txt.gz

**Sources**:
- From any Cell Ranger v9+ installation under `barcodes/` directory
- From 10x Genomics support documentation
- **Note**: Cell Ranger v9+ renamed the file with `_TRU` suffix, but both versions work

**What it is**: This is the set of known valid barcodes for the 3' v4 chemistry kit.

### 2. Reference Genome

You have two options for reference genomes:

#### Option A: Use 10x Prebuilt References (Recommended for Cell Ranger Parity)

**Download**: 10x Genomics → Support → Cell Ranger → Downloads → Reference packages

**Examples**:
- Human: `refdata-gex-GRCh38-2024-A.tar.gz`
- Mouse: `refdata-gex-GRCm39-2024-A.tar.gz`

**Advantages**:
- Highest parity with Cell Ranger results
- Pre-filtered to match 10x's gene annotation standards
- Includes optimized STAR indices

#### Option B: Build Custom References with 10x-like Filtering

Use WARP's `BuildIndices` pipeline with:
- GENCODE/Ensembl annotations
- 10x-style biotype filtering
- ID version stripping to match Cell Ranger format

## Usage

### Input Parameters

```json
{
  "Optimus.chemistry": "tenX_v4",
  "Optimus.whitelist": "gs://path/to/3M-3pgex-may-2023.txt.gz",
  "Optimus.tar_star_reference": "gs://path/to/star_reference.tar",
  "Optimus.annotations_gtf": "gs://path/to/annotations.gtf"
}
```

### Chemistry Validation

The pipeline automatically:
- Validates R1 read length (expects 28 bp for v4)
- Selects appropriate whitelist based on cloud provider
- Sets STARsolo parameters for 16-bp CB + 12-bp UMI

## Technical Details

### Read Structure
- **R1**: CB(1-16) + UMI(17-28) = 28 bp total
- **R2**: cDNA sequence

### STARsolo Parameters
```bash
--soloType CB_UMI_Simple
--soloCBstart 1 --soloCBlen 16
--soloUMIstart 17 --soloUMIlen 12
--soloCBwhitelist 3M-3pgex-may-2023.txt.gz
```

### Expected Results

With proper references, you should see:
- Cell overlap with Cell Ranger: ~73k+ cells (as in your comparison report)
- Per-cell Spearman correlation: ≥0.98
- Per-gene Spearman correlation: ≥0.97
- Reduced "present in Optimus, absent in CR" gene differences

## Troubleshooting

### Common Issues

1. **R1 length mismatch**: Ensure your FASTQ files have 28-bp R1 reads
2. **Whitelist not found**: Verify the `3M-3pgex-may-2023.txt.gz` file path
3. **Gene count differences**: Use 10x-style references for best parity

### Validation

Run on a small v4 dataset and verify:
- Chemistry detection selects "v4" branch
- STARsolo shows CB:1-16, UMI:17-28
- Whitelist resolves to 3M-3pgex-may-2023 file

## Example Files

See `example_inputs/human_v4_example.json` for a complete v4 configuration example.
