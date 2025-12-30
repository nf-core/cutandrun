# ChIPBinner: A Bin-Based Approach for Broad Histone Mark Analysis

## Overview

ChIPBinner takes a fundamentally different approach from traditional peak callers like MACS2, SEACR, or GoPeaks. Rather than defining discrete enriched regions, it divides the genome into uniform fixed-size bins (typically 10kb), quantifies signal in every bin, and uses unsupervised clustering (HDBSCAN) to identify bins with similar behavior across conditions.

**Key distinction:** No binary peak/non-peak classification—instead, it focuses on relative changes between conditions.

---

## When to Use [ChIPBinner](https://github.com/padilr1/ChIPbinner.git)

### Ideal Use Cases

- **Global changes** (e.g., NSD1 KO causing genome-wide H3K36me2 loss)—this is exactly what it was designed for
- **Broad histone marks** like H3K36me2, H3K36me3, H3K27me3
- When you have **external normalization** available (mass spec or spike-in)

### Less Suitable For

- **Sharp/punctate marks** (H3K4me3, H3K9ac, CTCF)—10kb bins are too large
- **Focal changes** affecting only a few loci
- When you need traditional peak calls as your primary output

### Assay Compatibility

Works with ChIP-seq, CUT&RUN, and CUT&Tag. For CUT&RUN/Tag, use IgG as your control instead of genomic input.

---

## Core Assumptions

1. **Most bins are NOT differentially enriched** (required for DESeq2/edgeR normalization)—explicitly violated with global changes, which reduces sensitivity

2. **Bin size captures relevant biology**—10kb works well for broad marks but misses fine structure in sharp marks

3. **Signal is comparable after normalization**—assumes your scaling factors (from MS or spike-in) are accurate

---

## Normalization Options

**Default normalization formula:**
```
log2( ((ChIP × scaling_factor) / (ChIP_depth/1e6) + pseudocount) /
      (Input / (Input_depth/1e6) + pseudocount) )
```

**Available methods:**

- **Library size** (`depth_norm=TRUE`)—CPM-like normalization
- **Input/control** (`use_input=TRUE`)—divides ChIP by input per bin; highly recommended
- **DESeq2 median ratio**—sample-specific size factors; assumes most bins unchanged
- **edgeR TMM**—trimmed mean of M-values; also assumes most bins unchanged
- **External scaling factors**—spike-in or mass spectrometry derived

---

## Handling Replicates

### Workflow for 2 Replicates × 2 Conditions

1. **Normalize each replicate individually** with its matched input via `norm_bw()`

2. **QC check:** Use `plot_PCA()` and `plot_correlation()` to verify replicate consistency

3. **Merge replicates per condition** using `merge_norm_bw()` for clustering and visualization

4. **Keep replicates separate for differential analysis**—ROTS requires individual replicates

```r
# Differential analysis needs individual replicates
differentialBinAnalysis(
  treated_sample_bigWigFiles = c("KO_rep1.bw", "KO_rep2.bw"),
  wildtype_sample_bigWigFiles = c("WT_rep1.bw", "WT_rep2.bw"),
  ...
)
```

---

## Differential Analysis with ROTS

ChIPBinner uses ROTS (Reproducibility-Optimized Test Statistic) rather than DESeq2 or edgeR.

**Why ROTS?** It outperforms other methods in datasets with a large proportion of differentially expressed features and skewed distributions—conditions frequently observed in ChIP-seq data following mutations affecting global histone levels.

**Output per bin:** logFC, FDR, and cluster assignment

---

## Unique Features

- **HDBSCAN clustering**—groups bins by behavior pattern (e.g., "bins that lose signal in KO")
- **Genic vs intergenic stratification**—separate analysis for gene bodies and intergenic regions
- **LOLA enrichment**—statistical overlap with functional annotations
- **Mass spectrometry normalization**—direct integration of MS-derived scaling factors
- **Exploratory tools**—PCA, correlation plots, density scatterplots

---

## Inputs & Outputs

### Required Inputs

- BAM files (aligned reads)
- Genome assembly (hg38 or mm10 only)
- Binning reference BED file

### Optional Inputs

- Control/Input BAM
- MS coefficients (YAML)
- Blacklist regions

### Key Outputs

- Normalized signal (bigWig)
- Signal matrix (CSV)
- HDBSCAN cluster assignments
- Annotated cluster GRanges objects
- Differential statistics per bin

---

## Snakemake Workflow Extension

The [snakemake_chipbinner](https://github.com/andygglez/snakemake_chipbinner) repo adds automation:

- Automated BAM → binned BED conversion
- CSV-based sample sheet management
- Config-driven MS normalization
- Batch comparisons
- Clustering parameter grid search (100–1000 for minpts/minsamps)

**Note:** Controls are disabled by default (`use_input: False`) in the Snakemake config—you may want to change this.

---
## [ChIPBinner Database](https://github.com/padilr1/ChIPbinner_database?tab=readme-ov-file)

The ChIPBinner database repo contains example data as well as curated databases from Ensembl, Encode and RepeatMasker for use with ChIPbinner.

'reference_windows' contains reference windows that can be used to transform aligned BAM files into binned BED files.

'example_data' contains complete datasets from Farhangdoost et al. (https://doi.org/10.1016/j.celrep.2021.108769):

HNSCC H3K36me2 ChIP-seq samples binned into 10kb windows
normalized bigWig files
complete pooled BED file of genomic coordinates
complete matrix file of binned scores
complete clustering results from running HDBSCAN on the whole dataset
'functional_db' contains curated databases (as R objects) for input into the 'enrich_clust()' function of ChIPbinner.

ensemblDB from https://useast.ensembl.org/info/genome/index.html"
ccreDB (candidate cis-regulatory elements) from "https://useast.ensembl.org/info/genome/index.html"
repeatsDB (repeatMasker) from https://www.repeatmasker.org/


---
## Summary

ChIPBinner is purpose-built for analyzing broad histone marks under conditions of global change. It's excellent for studying marks like H3K36me2/me3 after writer/eraser knockouts when you have external normalization available. For sharp marks or focal changes, stick with traditional peak callers.
