# Peak-Calling Extension Module Specification

## Overview

This specification defines a multi-caller peak-calling extension module for CUT&RUN data analysis. The module runs **7 peak-calling variants** across 3 tools (MACS2, epic2, SPAN/omnipeaks) on all samples, generates consensus peaks, computes QC metrics, and produces DiffBind-ready sample sheets.

---

## Table of Contents

1. [Peak-Calling Variants](#1-peak-calling-variants)
2. [Input Requirements](#2-input-requirements)
3. [Experimental Design](#3-experimental-design)
4. [Genome Handling](#4-genome-handling)
5. [Control/IgG Handling](#5-controligg-handling)
6. [Replicate Handling](#6-replicate-handling)
7. [Execution Model](#7-execution-model)
8. [Output Structure](#8-output-structure)
9. [QC Metrics and Reporting](#9-qc-metrics-and-reporting)
10. [DiffBind Integration](#10-diffbind-integration)
11. [Resource Requirements](#11-resource-requirements)
12. [Error Handling](#12-error-handling)
13. [Special Considerations](#13-special-considerations)

---

## 1. Peak-Calling Variants

### 1.1 MACS2 (2 variants) - No Control

CUT&RUN-optimized peak calling without IgG background subtraction.

| Variant | Command |
|---------|---------|
| `macs2_narrow` | `macs2 callpeak --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05 --gsize $GSIZE --format BAMPE --name $PREFIX --treatment $TREATMENT_BAM` |
| `macs2_broad` | `macs2 callpeak --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05 --broad --broad-cutoff 0.1 --gsize $GSIZE --format BAMPE --name $PREFIX --treatment $TREATMENT_BAM` |

**Key Parameters:**
- `--nomodel --shift -75 --extsize 150`: Fixed 150bp fragment handling (CUT&RUN optimized)
- `--keep-dup all`: Relies on upstream deduplication
- `-q 0.05`: FDR threshold (same for narrow and broad)
- `--broad-cutoff 0.1`: Linking threshold for broad peaks
- **No `--control`**: Peaks called without IgG background

**Environment:**
```bash
# Run in subshell to isolate Python 2.7 environment
(
  module purge
  module load viz
  module load biology
  module load py-macs2/2.1.1_py27
  module load samtools/1.8
  # ... run macs2 command ...
)
```

**Outputs:**
- Narrow: `{PREFIX}_peaks.narrowPeak`
- Broad: `{PREFIX}_peaks.broadPeak`
- Intermediates: `{PREFIX}_model.r`, `{PREFIX}_peaks.xls`, `{PREFIX}_summits.bed`

---

### 1.2 epic2 (3 variants) - With Control

SICER2-replacement for broad domain calling with IgG background subtraction.

| Variant | Window | Gaps | Command |
|---------|--------|------|---------|
| `epic2_200bp` | 200 | 3 (=600bp) | `epic2 --treatment $TREATMENT_BAM --control $POOLED_CONTROL_BAM --genome $GENOME -w 200 -g 3 -fdr 0.01 --output $PREFIX_epic2_200bp.peaks` |
| `epic2_150bp` | 150 | 2 (=300bp) | `epic2 --treatment $TREATMENT_BAM --control $POOLED_CONTROL_BAM --genome $GENOME -w 150 -g 2 -fdr 0.01 --output $PREFIX_epic2_150bp.peaks` |
| `epic2_25bp` | 25 | 2 (=50bp) | `epic2 --treatment $TREATMENT_BAM --control $POOLED_CONTROL_BAM --genome $GENOME -w 25 -g 2 -fdr 0.01 --output $PREFIX_epic2_25bp.peaks` |

**Key Parameters:**
- `-w/--window-size`: Bin size in base pairs
- `-g/--gaps-allowed`: Maximum gaps between enriched windows (measured in windows, not bp)
- `-fdr 0.01`: FDR threshold
- `--control`: Uses pooled IgG control (see Section 5)

**Environment:**
```bash
(
  module purge
  module load python/3.12.1
  # ... run epic2 command ...
)
```

**Outputs:**
- Main: `{PREFIX}_epic2_{variant}.peaks`

---

### 1.3 SPAN/omnipeaks (2 variants) - With Control

Semi-parametric peak calling, particularly suited for broad histone marks.

| Variant | FDR | Gap | Command |
|---------|-----|-----|---------|
| `span_default` | 0.05 | 5 | `java -Xmx8G -jar $OMNIPEAKS_JAR analyze -t $TREATMENT_BAM -c $POOLED_CONTROL_BAM --cs $CHROM_SIZES -p $PREFIX_span_default` |
| `span_stringent` | 1.0E-9 | 2 | `java -Xmx8G -jar $OMNIPEAKS_JAR analyze -t $TREATMENT_BAM -c $POOLED_CONTROL_BAM --cs $CHROM_SIZES --fdr 1.0E-9 --gap 2 -p $PREFIX_span_stringent` |

**Key Parameters:**
- `--bin`: Window size (default 200bp)
- `--fdr`: FDR threshold
- `--gap`: Maximum gaps between enriched bins
- `--cs`: Chromosome sizes file (generated from BAM, see Section 4.3)
- `-Xmx8G`: 8GB Java heap memory

**Environment:**
```bash
(
  module purge
  module load java/21.0.4
  export OMNIPEAKS_JAR=/home/groups/ogozani/programs/omnipeak/omnipeak.jar
  # ... run SPAN command ...
)
```

**Outputs:**
- Peak file: `{PREFIX}.peak`

---

## 2. Input Requirements

### 2.1 BAM Files

All commands expect nf-core/cutandrun style BAM files:
- **Treatment:** `.target.markdup.sorted.bam` from alignment output
- **Control/IgG:** `.target.dedup.sorted.bam` from alignment output

### 2.2 BAM Index Validation

Before peak calling, validate that BAM index files exist:
```bash
for bam in *.bam; do
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    samtools index "$bam"
  fi
done
```

### 2.3 Sample Sheet

Use the nf-core samplesheet.csv format. Reference: [nf-core/cutandrun samplesheet documentation](https://nf-co.re/cutandrun/usage#samplesheet-input)

**Required columns:**
- `sample`: Sample identifier
- `group`: Format `{mark}_{condition}` (e.g., `H3K36me2_WT`, `H3K4me3_KO`)
- `replicate`: Replicate number (1 or 2)
- `fastq_1`, `fastq_2`: FASTQ paths (for upstream pipeline)
- `control`: IgG control sample name

**Validation:**
- Parse `group` column to extract mark and condition
- Validate against BAM filename patterns
- Ensure 2×2 structure per mark (see Section 3)

---

## 3. Experimental Design

### 3.1 Enforced Structure

The pipeline enforces a **2×2 design** per histone mark:
- **2 conditions**: e.g., WT (wild-type/control) and KO (knockout/treated)
- **2 biological replicates** per condition

```
Per mark:
├── Condition_A (e.g., WT)
│   ├── Replicate 1
│   └── Replicate 2
└── Condition_B (e.g., KO)
    ├── Replicate 1
    └── Replicate 2
```

### 3.2 Validation

At pipeline start, validate that each mark has exactly:
- 2 unique conditions
- 2 replicates per condition
- 4 total samples per mark

Fail with descriptive error if structure is not met.

---

## 4. Genome Handling

### 4.1 Genome Detection

Auto-detect genome build from BAM headers with user validation:

```bash
# Extract genome info from BAM header
samtools view -H "$BAM" | grep "^@SQ" | head -1

# Check for UCSC notation (chr prefix)
if samtools view -H "$BAM" | grep -q "SN:chr"; then
  CHR_NOTATION="UCSC"
else
  CHR_NOTATION="Ensembl"
fi

# Detect species by contig lengths
# chr1/1 length ~248M = human (hg38)
# chr1/1 length ~195M = mouse (mm39/mm10)
```

If ambiguous, prompt user to confirm genome build.

### 4.2 Chromosome Notation Conversion

If BAM files use Ensembl notation (1, 2, ...) but tools require UCSC (chr1, chr2, ...), auto-convert:

```bash
# Convert Ensembl to UCSC notation
samtools view -H "$BAM" | \
  sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' \
      -e 's/SN:MT/SN:chrM/' | \
  samtools reheader - "$BAM" > "${BAM%.bam}_ucsc.bam"
samtools index "${BAM%.bam}_ucsc.bam"
```

**Storage:** Save converted BAMs permanently alongside originals for reuse.

### 4.3 Genome Size Configuration

| Species | Genome | MACS2 `--gsize` | epic2/SPAN `--genome` |
|---------|--------|-----------------|----------------------|
| Human | hg38/GRCh38 | `2.7E+9` | `hg38` |
| Mouse | mm39/GRCm39 | `1.87E+9` | `mm39` → `mm10`* |

*Note: Auto-map mm39 to mm10 for tools that don't recognize mm39 (similar enough for peak calling).

### 4.4 Chromosome Sizes Generation

Generate chrom.sizes from BAM header for SPAN:

```bash
samtools view -H "$BAM" | \
  grep "^@SQ" | \
  sed 's/@SQ\tSN:\([^\t]*\)\tLN:\([0-9]*\).*/\1\t\2/' > "${GENOME}.chrom.sizes"
```

---

## 5. Control/IgG Handling

### 5.1 Control Usage by Tool

| Tool | Uses Control | Control Source |
|------|--------------|----------------|
| MACS2 | **No** | N/A |
| epic2 | **Yes** | Pooled IgG |
| SPAN | **Yes** | Pooled IgG |

### 5.2 Pooled IgG Creation

**Strategy:** Pre-create pooled IgG BAM once at pipeline start, reuse for all callers.

**Pooling method:** Simple merge with samtools (no depth normalization):

```bash
# Pool IgG controls across replicates within each condition
samtools merge -@ 8 \
  "${CONDITION}_pooled_IgG.bam" \
  "${CONDITION}_rep1_IgG.bam" \
  "${CONDITION}_rep2_IgG.bam"
samtools index "${CONDITION}_pooled_IgG.bam"
```

### 5.3 Missing Control Fallback

If a sample's matched IgG fails QC or is missing:
1. Log warning: "IgG control for {sample} unavailable, using pooled IgG from {other_condition}"
2. Use pooled IgG from the other experimental condition
3. Flag in QC report

---

## 6. Replicate Handling

### 6.1 Peak Calling Scope

Peak calling runs on **individual replicates only** (not pooled).

### 6.2 Consensus Peak Generation

After per-replicate peak calling, generate consensus peaks for each caller:

**Intersection (conservative):** Peaks present in both replicates with ≥20% overlap:
```bash
bedtools intersect -a rep1_peaks.bed -b rep2_peaks.bed -f 0.2 -wa > consensus_peaks.bed
```

**Union (permissive):** All peaks from either replicate, merged:
```bash
cat rep1_peaks.bed rep2_peaks.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge > union_peaks.bed
```

### 6.3 Replicate Imbalance Warning

If one replicate has >3× more peaks than the other, log warning:
```
WARNING: Replicate imbalance detected for {sample}
  Rep1: 15,234 peaks
  Rep2: 4,102 peaks (3.7x difference)
  Consider investigating replicate quality.
```

---

## 7. Execution Model

### 7.1 Parallel Execution

- **All 7 variants** run on every sample
- **Concurrent callers** allowed on same sample
- Each caller runs in its own **subshell** with isolated module environment

### 7.2 Module Isolation

To handle Python version conflicts (MACS2 requires Py2.7, epic2 requires Py3.12):

```bash
# Each caller runs in isolated subshell
run_macs2() {
  (
    module purge
    module load py-macs2/2.1.1_py27
    macs2 callpeak ...
  )
}

run_epic2() {
  (
    module purge
    module load python/3.12.1
    epic2 ...
  )
}
```

### 7.3 Module Loading

**Primary approach:** Source existing modules.env file:
```bash
if [[ -f "/path/to/modules.env" ]]; then
  source "/path/to/modules.env"
else
  # Fallback: inline module loads
  module load R/4.4.2
  module load python/3.12.1
  module load biology
  # ... etc
fi
```

### 7.4 Dry-Run Mode

**Required feature:** Validate inputs and show planned execution without running:

```bash
./peak_calling.sh --dry-run --samplesheet samples.csv

# Output:
# [DRY-RUN] Detected genome: hg38 (UCSC notation)
# [DRY-RUN] Samples: 8 (2 conditions × 2 replicates × 2 marks)
# [DRY-RUN] Will create pooled IgG: WT_pooled_IgG.bam, KO_pooled_IgG.bam
# [DRY-RUN] Will run 7 peak-callers on each of 8 samples (56 total jobs)
# [DRY-RUN] Estimated resource usage: ...
```

---

## 8. Output Structure

### 8.1 Directory Layout

```
04_called_peaks/
├── macs2_narrow/
│   ├── {condition}_{mark}_{rep}_peaks.narrowPeak
│   └── ...
├── macs2_broad/
│   ├── {condition}_{mark}_{rep}_peaks.broadPeak
│   └── ...
├── epic2_200bp/
│   ├── {condition}_{mark}_{rep}_epic2_200bp.peaks
│   └── ...
├── epic2_150bp/
├── epic2_25bp/
├── span_default/
│   ├── {condition}_{mark}_{rep}_span_default.peak
│   └── ...
├── span_stringent/
├── consensus/
│   ├── macs2_narrow/
│   │   ├── {mark}_{condition}_consensus.bed
│   │   └── {mark}_{condition}_union.bed
│   └── ... (for each caller)
├── pooled_controls/
│   ├── {condition}_pooled_IgG.bam
│   └── {condition}_pooled_IgG.bam.bai
└── ucsc_converted/
    └── ... (UCSC-notation BAMs if converted)

05_intermediate_files/
├── macs2/
│   ├── {sample}_model.r
│   ├── {sample}_peaks.xls
│   └── ...
└── logs/
    ├── macs2_narrow_{sample}.log
    └── ... (full stdout/stderr per caller per sample)

06_reports/
├── peak_calling_summary.html
├── diffbind_samplesheets/
│   ├── macs2_narrow_samplesheet.csv
│   ├── macs2_broad_samplesheet.csv
│   ├── epic2_200bp_samplesheet.csv
│   └── ... (7 total)
└── qc_metrics/
    ├── frip_summary.tsv
    └── peak_counts.tsv
```

### 8.2 File Naming Convention

Format: `{condition}_{mark}_{rep}_{caller}_{suffix}`

Examples:
- `WT_H3K36me2_rep1_macs2_narrow_peaks.narrowPeak`
- `KO_H3K4me3_rep2_epic2_200bp.peaks`
- `WT_H3K36me2_rep1_span_stringent.peak`

---

## 9. QC Metrics and Reporting

### 9.1 FRiP Calculation

Calculate Fraction of Reads in Peaks for every sample × caller:

```bash
# Total reads in BAM
total_reads=$(samtools view -c -F 4 "$BAM")

# Reads overlapping peaks
reads_in_peaks=$(bedtools intersect -a "$BAM" -b "$PEAKS" -bed -u | wc -l)

# FRiP
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
```

### 9.2 Best Caller Recommendation

**Metric:** Balanced score = FRiP × log10(peak_count)

**Scope:** Per-mark across all samples (not per-sample)

**Threshold:** Flag callers with < 500 peaks as potentially failed.

### 9.3 HTML Report

Generate interactive HTML report using plotly with:

1. **FRiP Barplots:** Compare FRiP across all 7 callers, grouped by sample
2. **Peak Count Tables:** Tabular summary of peak counts per caller × sample
3. **Venn Diagrams:** Peak overlap between callers (for each sample)
4. **Heatmaps:** Pairwise peak overlap similarity between all caller pairs
5. **Best Caller Recommendation:** Highlighted per-mark recommendation

### 9.4 Logging

**Verbosity:** Capture full stdout/stderr from each tool.

Log files: `05_intermediate_files/logs/{caller}_{sample}.log`

---

## 10. DiffBind Integration

### 10.1 Sample Sheet Generation

Generate **separate DiffBind-ready sample sheets** for each of the 7 peak-callers.

**Format:** (matches DiffBind requirements)
```csv
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
WT_H3K36me2_rep1,CellLine,H3K36me2,WT,1,/path/to/WT_H3K36me2_rep1.bam,WT_IgG,/path/to/WT_pooled_IgG.bam,/path/to/peaks.narrowPeak,narrow
```

### 10.2 Peak File References

DiffBind sheets point to **per-replicate peaks** (not consensus).

DiffBind handles replicate merging internally during differential analysis.

### 10.3 Output Files

```
06_reports/diffbind_samplesheets/
├── macs2_narrow_samplesheet.csv
├── macs2_broad_samplesheet.csv
├── epic2_200bp_samplesheet.csv
├── epic2_150bp_samplesheet.csv
├── epic2_25bp_samplesheet.csv
├── span_default_samplesheet.csv
└── span_stringent_samplesheet.csv
```

---

## 11. Resource Requirements

### 11.1 Per-Caller Resources

| Caller | Memory | CPUs | Time (est.) |
|--------|--------|------|-------------|
| macs2_narrow | 4 GB | 1 | 5-15 min |
| macs2_broad | 4 GB | 1 | 5-15 min |
| epic2_200bp | 8 GB | 8 | 10-30 min |
| epic2_150bp | 8 GB | 8 | 10-30 min |
| epic2_25bp | 16 GB | 8 | 20-60 min |
| span_default | 8 GB | 1 | 15-45 min |
| span_stringent | 8 GB | 1 | 15-45 min |

### 11.2 SLURM Job Template

```bash
#SBATCH --job-name=peakcall_{caller}
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem={MEM}
#SBATCH --cpus-per-task={CPUS}
```

---

## 12. Error Handling

### 12.1 Failure Strategy

**Continue on failure:** If a peak-caller fails on a sample:
1. Log the failure with exit code and stderr
2. Continue with other callers on that sample
3. Continue with other samples
4. Report all failures in final summary

### 12.2 Failure Logging

```
ERROR: macs2_narrow failed on WT_H3K36me2_rep1
  Exit code: 1
  Stderr: [captured stderr]
  Log file: 05_intermediate_files/logs/macs2_narrow_WT_H3K36me2_rep1.log
```

### 12.3 Summary Report

At pipeline end, generate failure summary:
```
Peak Calling Summary
====================
Total jobs: 56 (8 samples × 7 callers)
Successful: 54
Failed: 2
  - macs2_narrow on WT_H3K36me2_rep1: exit 1
  - span_stringent on KO_H3K4me3_rep2: exit 137 (OOM)
```

---

## 13. Special Considerations

### 13.1 Spike-In Normalized Data

When working with spike-in normalized BAMs:
- **MACS2:** Consider using `--scale-to small` or adjusting `--extsize`
- **epic2/SPAN:** Spike-in normalization is typically applied upstream; no parameter changes needed
- **FRiP calculation:** Use spike-in normalized read counts if available

### 13.2 MS-Normalized Data

For mass-spectrometry normalized data:
- Normalization factors are typically applied to coverage/signal tracks, not peak calling
- Peak calling can proceed on raw BAMs
- Document normalization factors in metadata for downstream analysis

### 13.3 Low-Depth Samples

No minimum depth requirement enforced, but:
- Peak callers may produce few or no peaks on very low-depth samples
- FRiP will naturally reflect low signal
- The < 500 peaks warning will flag potentially problematic samples

---

## Appendix A: Complete Command Reference

### A.1 MACS2 Narrow
```bash
macs2 callpeak \
  --nomodel \
  --shift -75 \
  --extsize 150 \
  --keep-dup all \
  -q 0.05 \
  --gsize 2.7E+9 \
  --format BAMPE \
  --name ${PREFIX} \
  --treatment ${TREATMENT_BAM} \
  --outdir ${OUTDIR}
```

### A.2 MACS2 Broad
```bash
macs2 callpeak \
  --nomodel \
  --shift -75 \
  --extsize 150 \
  --keep-dup all \
  -q 0.05 \
  --broad \
  --broad-cutoff 0.1 \
  --gsize 2.7E+9 \
  --format BAMPE \
  --name ${PREFIX} \
  --treatment ${TREATMENT_BAM} \
  --outdir ${OUTDIR}
```

### A.3 epic2 (200bp example)
```bash
epic2 \
  --treatment ${TREATMENT_BAM} \
  --control ${POOLED_CONTROL_BAM} \
  --genome hg38 \
  --window-size 200 \
  --gaps-allowed 3 \
  --false-discovery-rate-cutoff 0.01 \
  --output ${PREFIX}_epic2_200bp.peaks
```

### A.4 SPAN Default
```bash
java -Xmx8G -jar ${OMNIPEAKS_JAR} analyze \
  -t ${TREATMENT_BAM} \
  -c ${POOLED_CONTROL_BAM} \
  --cs ${CHROM_SIZES} \
  -p ${PREFIX}_span_default
```

### A.5 SPAN Stringent
```bash
java -Xmx8G -jar ${OMNIPEAKS_JAR} analyze \
  -t ${TREATMENT_BAM} \
  -c ${POOLED_CONTROL_BAM} \
  --cs ${CHROM_SIZES} \
  --fdr 1.0E-9 \
  --gap 2 \
  -p ${PREFIX}_span_stringent
```

---

## Appendix B: Environment Reference

### B.1 modules.env Location
```
/scratch/groups/ogozani/seq_data/dylan/sirt7/maria/23_01_13_cnr_sirt7/nf-core_output_sirt7_A549_230113_bpmNorm/clones/worktree_refactor/05_results_00_AddNewPeakCallers/agents/modules.env
```

### B.2 Tool Locations
```bash
# SPAN/omnipeaks
export OMNIPEAKS_JAR=/home/groups/ogozani/programs/omnipeak/omnipeak.jar

# ChromHMM (if needed for downstream)
export CHROMHMM_JAR=/home/groups/ogozani/programs/ChromHMM/ChromHMM.jar
```

### B.3 Required Modules
```bash
module load R/4.4.2           # For downstream analysis
module load python/3.12.1     # For epic2, general Python
module load biology           # Bioinformatics tools
module load py-deeptools/3.5.6_py312
module load bedtools/2.30.0   # Consensus peaks
module load samtools/1.16.1   # BAM operations
module load parallel/20200822 # Parallel execution
module load bedops/2.4.40     # Fast set operations
module load java/21.0.4       # For SPAN
```

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2024-12-29 | Claude Code | Initial specification |
