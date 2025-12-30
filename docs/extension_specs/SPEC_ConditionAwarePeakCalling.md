# SPEC: Condition-Aware Peak Calling Extension

## Overview

This specification defines extensions to nf-core/cutandrun that add:
1. **Condition metadata** in the samplesheet for experimental design (Treatment vs. Control)
2. **Group-specific handling** for normalization and peak calling
3. **Support for additional peak callers**: GoPeaks, epic2, and SPAN/omnipeaks

These features enable differential analysis workflows and provide users with flexibility in choosing peak-calling strategies optimized for different histone mark profiles.

---

## Table of Contents

1. [Motivation](#1-motivation)
2. [Samplesheet Extension](#2-samplesheet-extension)
3. [Group-Specific Processing](#3-group-specific-processing)
4. [Additional Peak Callers](#4-additional-peak-callers)
5. [Implementation Plan](#5-implementation-plan)
6. [Configuration Parameters](#6-configuration-parameters)
7. [Output Structure](#7-output-structure)
8. [Testing Strategy](#8-testing-strategy)
9. [Migration Guide](#9-migration-guide)

---

## 1. Motivation

### 1.1 Current Limitations

The current nf-core/cutandrun pipeline:
- Treats all non-control samples as equivalent targets
- Applies uniform normalization across all samples
- Supports only SEACR and MACS2 for peak calling
- Cannot distinguish between experimental conditions (e.g., WT vs. KO)

### 1.2 Use Cases Enabled

This extension supports:
- **Differential analysis preparation**: Generate DiffBind-ready samplesheets with condition metadata
- **Condition-aware normalization**: Apply appropriate scaling factors per experimental group
- **Peak caller comparison**: Run multiple peak callers to identify optimal methods for specific marks
- **Broad histone mark analysis**: epic2 and SPAN are better suited for broad domains (H3K36me2/me3, H3K27me3)

---

## 2. Samplesheet Extension

### 2.1 New Samplesheet Format

**Current format:**
```csv
group,replicate,fastq_1,fastq_2,control
```

**Extended format:**
```csv
group,replicate,fastq_1,fastq_2,control,condition
```

### 2.2 Column Definitions

| Column | Required | Description | Example Values |
|--------|----------|-------------|----------------|
| `group` | Yes | Sample group identifier (typically histone mark) | `H3K36me2`, `H3K4me3` |
| `replicate` | Yes | Replicate number (integer ≥ 1) | `1`, `2` |
| `fastq_1` | Yes | Path to read 1 FASTQ file | `/path/to/sample_R1.fq.gz` |
| `fastq_2` | Yes | Path to read 2 FASTQ file | `/path/to/sample_R2.fq.gz` |
| `control` | Conditional | IgG control group name (empty for controls) | `IgG_WT`, `IgG_KO` |
| `condition` | No* | Experimental condition | `WT`, `KO`, `Treatment`, `Control` |

*`condition` is optional for backward compatibility but required for condition-aware processing.

### 2.3 Example Samplesheet

```csv
group,replicate,fastq_1,fastq_2,control,condition
H3K36me2,1,H3K36me2_WT_rep1_R1.fq.gz,H3K36me2_WT_rep1_R2.fq.gz,IgG_WT,WT
H3K36me2,2,H3K36me2_WT_rep2_R1.fq.gz,H3K36me2_WT_rep2_R2.fq.gz,IgG_WT,WT
H3K36me2,1,H3K36me2_KO_rep1_R1.fq.gz,H3K36me2_KO_rep1_R2.fq.gz,IgG_KO,KO
H3K36me2,2,H3K36me2_KO_rep2_R1.fq.gz,H3K36me2_KO_rep2_R2.fq.gz,IgG_KO,KO
H3K4me3,1,H3K4me3_WT_rep1_R1.fq.gz,H3K4me3_WT_rep1_R2.fq.gz,IgG_WT,WT
H3K4me3,2,H3K4me3_WT_rep2_R1.fq.gz,H3K4me3_WT_rep2_R2.fq.gz,IgG_WT,WT
H3K4me3,1,H3K4me3_KO_rep1_R1.fq.gz,H3K4me3_KO_rep1_R2.fq.gz,IgG_KO,KO
H3K4me3,2,H3K4me3_KO_rep2_R1.fq.gz,H3K4me3_KO_rep2_R2.fq.gz,IgG_KO,KO
IgG_WT,1,IgG_WT_rep1_R1.fq.gz,IgG_WT_rep1_R2.fq.gz,,WT
IgG_WT,2,IgG_WT_rep2_R1.fq.gz,IgG_WT_rep2_R2.fq.gz,,WT
IgG_KO,1,IgG_KO_rep1_R1.fq.gz,IgG_KO_rep1_R2.fq.gz,,KO
IgG_KO,2,IgG_KO_rep2_R1.fq.gz,IgG_KO_rep2_R2.fq.gz,,KO
```

### 2.4 Validation Rules

The extended `check_samplesheet.py` must enforce:

1. **Condition consistency**: All replicates of a `group` × `condition` combination must have matching control assignments
2. **Condition format**: No spaces allowed; alphanumeric with underscores permitted
3. **Backward compatibility**: If `condition` column is absent or empty, default to `"default"` condition
4. **Control conditions**: IgG controls should have a condition that matches their target samples

### 2.5 Metadata Propagation

The `condition` field propagates through the pipeline via the `meta` map:

```groovy
meta.id            = row.id
meta.group         = row.group
meta.replicate     = row.replicate.toInteger()
meta.single_end    = row.single_end.toBoolean()
meta.is_control    = row.is_control.toBoolean()
meta.control_group = meta.is_control ? meta.group : row.control
meta.condition     = row.condition ?: "default"  // NEW FIELD
```

---

## 3. Group-Specific Processing

### 3.1 Condition-Aware Normalization

#### 3.1.1 Current Behavior

Currently, `PREPARE_PEAKCALLING` applies uniform normalization:
- Spike-in normalization: scale factor = `normalisation_c / spikein_aligned_reads`
- Other modes: RPKM, CPM, BPM, RPGC, or None

#### 3.1.2 Extended Behavior

Add a new parameter `--condition_aware_norm` (default: `false`) that enables:

**Per-condition spike-in normalization:**
```groovy
// Calculate condition-specific scaling
ch_metadata.branch { it ->
    condition_A: it.condition == "WT"
    condition_B: it.condition == "KO"
}
.set { ch_metadata_by_condition }

// Compute median spike-in count per condition for more robust normalization
// Scale each sample relative to its condition's median
```

**Use case:** When global histone levels change between conditions (e.g., NSD1 KO reducing H3K36me2), condition-specific normalization prevents artifacts.

#### 3.1.3 Normalization Parameter Extension

New parameter: `--normalisation_scope`

| Value | Behavior |
|-------|----------|
| `global` | Current behavior: all samples normalized together (default) |
| `condition` | Samples normalized within their condition group |
| `group` | Samples normalized within their `group` (histone mark) |
| `group_condition` | Samples normalized within `group` × `condition` |

### 3.2 Condition-Aware Peak Calling

#### 3.2.1 Control Pairing

Extend control pairing logic to match by condition:

**Current logic:**
```groovy
ch_bam_control.map{ row -> [row[0].control_group + "_" + row[0].replicate, row] }
.cross( ch_bam_target.map{ row -> [row[0].control_group, row] } )
```

**Extended logic:**
```groovy
// Match controls to targets by both control_group and condition
ch_bam_control.map{ row ->
    [row[0].control_group + "_" + row[0].replicate + "_" + row[0].condition, row]
}
.cross( ch_bam_target.map{ row ->
    [row[0].control_group + "_" + row[0].condition, row]
} )
```

#### 3.2.2 IgG Pooling Strategy

For peak callers requiring IgG controls (epic2, SPAN), pool IgG **per condition**:

```groovy
// Pool IgG controls by condition
ch_bam_control
    .map { row -> [row[0].condition, row[1]] }
    .groupTuple()
    .map { condition, bams ->
        def new_meta = [id: "${condition}_pooled_IgG", condition: condition, is_control: true]
        [new_meta, bams.flatten()]
    }
    .set { ch_pooled_igg_by_condition }

// Merge BAMs for pooled control
SAMTOOLS_MERGE(ch_pooled_igg_by_condition)
```

#### 3.2.3 Consensus Peaks by Condition

Extend `consensus_peak_mode` parameter:

| Value | Behavior |
|-------|----------|
| `group` | Consensus within `group` (current default) |
| `all` | Single consensus from all samples (current) |
| `condition` | Consensus within each `condition` |
| `group_condition` | Consensus within `group` × `condition` combinations |

---

## 4. Additional Peak Callers

### 4.1 Supported Peak Callers

| Caller | Variants | Control Required | Best For |
|--------|----------|------------------|----------|
| SEACR | (existing) | Optional | Punctate marks |
| MACS2 | narrow, broad | Optional | General purpose |
| GoPeaks | narrow, broad | Optional | CUT&RUN optimized |
| epic2 | 200bp, 150bp, 25bp | Yes (pooled) | Broad domains |
| SPAN | default, stringent | Yes (pooled) | Broad histone marks |

### 4.2 Parameter: `--peakcaller`

Extend to accept comma-separated list of callers:

```bash
--peakcaller 'seacr,macs2_narrow,macs2_broad,gopeaks_narrow,gopeaks_broad,epic2_200bp,epic2_150bp,epic2_25bp,span_default,span_stringent'
```

**Validation:**
```groovy
def valid_callers = [
    'seacr',
    'macs2', 'macs2_narrow', 'macs2_broad',
    'gopeaks', 'gopeaks_narrow', 'gopeaks_broad',
    'epic2', 'epic2_200bp', 'epic2_150bp', 'epic2_25bp',
    'span', 'span_default', 'span_stringent'
]

callers = params.peakcaller.split(',').collect{ it.trim().toLowerCase() }
def invalid = callers - valid_callers
if (invalid) {
    exit 1, "Invalid peak caller(s): ${invalid.join(', ')}. Valid options: ${valid_callers.join(', ')}"
}
```

### 4.3 GoPeaks Module

#### 4.3.1 Module Definition

**File:** `modules/nf-core/gopeaks/main.nf`

```groovy
process GOPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gopeaks=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gopeaks:1.0.0--h9ee0642_0' :
        'biocontainers/gopeaks:1.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(control_bam)

    output:
    tuple val(meta), path("*_peaks.bed"), emit: peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam ? "-c $control_bam" : ''
    """
    gopeaks \\
        -b $bam \\
        $control \\
        -o ${prefix}_peaks.bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gopeaks: \$(gopeaks --version 2>&1 | sed 's/gopeaks version //')
    END_VERSIONS
    """
}
```

#### 4.3.2 GoPeaks Parameters

| Parameter | Narrow | Broad |
|-----------|--------|-------|
| `--mdist` | 1000 (default) | 3000 |
| `--pval` | 0.05 | 0.05 |
| `--step` | 100 | 500 |
| `--slide` | 50 | 250 |

**Module config (`conf/modules.config`):**
```groovy
withName: 'GOPEAKS_NARROW' {
    ext.args = '--mdist 1000 --pval 0.05 --step 100 --slide 50'
    ext.prefix = { "${meta.id}_gopeaks_narrow" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/gopeaks_narrow" },
        mode: params.publish_dir_mode
    ]
}

withName: 'GOPEAKS_BROAD' {
    ext.args = '--mdist 3000 --pval 0.05 --step 500 --slide 250'
    ext.prefix = { "${meta.id}_gopeaks_broad" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/gopeaks_broad" },
        mode: params.publish_dir_mode
    ]
}
```

### 4.4 epic2 Module

#### 4.4.1 Module Definition

**File:** `modules/nf-core/epic2/main.nf`

```groovy
process EPIC2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::epic2=0.0.52"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epic2:0.0.52--py39h14c64f4_0' :
        'biocontainers/epic2:0.0.52--py39h14c64f4_0' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    val genome

    output:
    tuple val(meta), path("*.peaks"), emit: peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epic2 \\
        --treatment $treatment_bam \\
        --control $control_bam \\
        --genome $genome \\
        --output ${prefix}.peaks \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epic2: \$(epic2 --version 2>&1 | sed 's/epic2 //')
    END_VERSIONS
    """
}
```

#### 4.4.2 epic2 Variants

| Variant | Window (`-w`) | Gaps (`-g`) | Gap (bp) | Use Case |
|---------|---------------|-------------|----------|----------|
| `epic2_200bp` | 200 | 3 | 600 | Default broad |
| `epic2_150bp` | 150 | 2 | 300 | Moderate resolution |
| `epic2_25bp` | 25 | 2 | 50 | High resolution |

**Module config:**
```groovy
withName: 'EPIC2_200BP' {
    ext.args = '--window-size 200 --gaps-allowed 3 --false-discovery-rate-cutoff 0.01'
    ext.prefix = { "${meta.id}_epic2_200bp" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/epic2_200bp" },
        mode: params.publish_dir_mode
    ]
}

withName: 'EPIC2_150BP' {
    ext.args = '--window-size 150 --gaps-allowed 2 --false-discovery-rate-cutoff 0.01'
    ext.prefix = { "${meta.id}_epic2_150bp" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/epic2_150bp" },
        mode: params.publish_dir_mode
    ]
}

withName: 'EPIC2_25BP' {
    ext.args = '--window-size 25 --gaps-allowed 2 --false-discovery-rate-cutoff 0.01'
    ext.prefix = { "${meta.id}_epic2_25bp" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/epic2_25bp" },
        mode: params.publish_dir_mode
    ]
}
```

### 4.5 SPAN/omnipeaks Module

#### 4.5.1 Module Definition

**File:** `modules/local/span/main.nf`

```groovy
process SPAN {
    tag "$meta.id"
    label 'process_high_memory'

    container 'docker://ghcr.io/jetbrains-research/span:latest'

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    path chrom_sizes

    output:
    tuple val(meta), path("*.peak"), emit: peaks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam.name != 'NO_CONTROL' ? "-c $control_bam" : ''
    def memory = task.memory.toGiga()
    """
    java -Xmx${memory}G -jar /opt/span/span.jar analyze \\
        -t $treatment_bam \\
        $control \\
        --cs $chrom_sizes \\
        -p ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        span: \$(java -jar /opt/span/span.jar --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
```

#### 4.5.2 SPAN Variants

| Variant | FDR | Gap | Use Case |
|---------|-----|-----|----------|
| `span_default` | 0.05 | 5 | General broad marks |
| `span_stringent` | 1.0E-9 | 2 | High-confidence peaks |

**Module config:**
```groovy
withName: 'SPAN_DEFAULT' {
    ext.args = '--fdr 0.05 --gap 5 --bin 200'
    ext.prefix = { "${meta.id}_span_default" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/span_default" },
        mode: params.publish_dir_mode
    ]
}

withName: 'SPAN_STRINGENT' {
    ext.args = '--fdr 1.0E-9 --gap 2 --bin 200'
    ext.prefix = { "${meta.id}_span_stringent" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/span_stringent" },
        mode: params.publish_dir_mode
    ]
}
```

### 4.6 MACS2 Extended Support

#### 4.6.1 Narrow and Broad Variants

Extend existing MACS2 module to explicitly support both modes:

| Variant | Key Parameters |
|---------|----------------|
| `macs2_narrow` | `--nomodel --shift -75 --extsize 150 -q 0.05` |
| `macs2_broad` | `--nomodel --shift -75 --extsize 150 -q 0.05 --broad --broad-cutoff 0.1` |

**Note:** Parameters are CUT&RUN-optimized (150bp fragment handling, no model building).

**Module config:**
```groovy
withName: 'MACS2_CALLPEAK_NARROW' {
    ext.args = '--nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05'
    ext.prefix = { "${meta.id}_macs2_narrow" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/macs2_narrow" },
        mode: params.publish_dir_mode
    ]
}

withName: 'MACS2_CALLPEAK_BROAD' {
    ext.args = '--nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05 --broad --broad-cutoff 0.1'
    ext.prefix = { "${meta.id}_macs2_broad" }
    publishDir = [
        path: { "${params.outdir}/04_called_peaks/macs2_broad" },
        mode: params.publish_dir_mode
    ]
}
```

---

## 5. Implementation Plan

### 5.1 Phase 1: Samplesheet Extension

**Files to modify:**

1. **`bin/check_samplesheet.py`**
   - Add `condition` column parsing
   - Implement backward compatibility (default condition)
   - Add validation for condition consistency
   - Update output header to include `condition`

2. **`subworkflows/local/input_check.nf`**
   - Update `get_samplesheet_paths()` to extract `condition`
   - Add `meta.condition` to metadata map

3. **`nextflow_schema.json`**
   - Document new samplesheet format
   - Add `condition` column description

### 5.2 Phase 2: Group-Specific Processing

**Files to modify:**

1. **`nextflow.config`**
   - Add `normalisation_scope` parameter
   - Add `condition_aware_norm` parameter
   - Add `consensus_peak_mode` options

2. **`subworkflows/local/prepare_peakcalling.nf`**
   - Implement condition-aware spike-in normalization
   - Add branching logic for different normalization scopes

3. **`subworkflows/local/consensus_peaks.nf`**
   - Extend grouping logic for condition-based consensus

4. **`workflows/cutandrun.nf`**
   - Update control pairing to respect conditions
   - Implement pooled IgG creation per condition

### 5.3 Phase 3: Additional Peak Callers

**New files:**

1. **`modules/nf-core/gopeaks/main.nf`** - GoPeaks process
2. **`modules/nf-core/epic2/main.nf`** - epic2 process
3. **`modules/local/span/main.nf`** - SPAN process
4. **`subworkflows/local/pool_controls.nf`** - IgG pooling subworkflow

**Files to modify:**

1. **`workflows/cutandrun.nf`**
   - Add includes for new modules
   - Implement multi-caller orchestration
   - Update primary/secondary peak assignment logic

2. **`conf/modules.config`**
   - Add configuration for all new caller variants

3. **`conf/base.config`**
   - Add resource labels for epic2/SPAN (high memory)

4. **`nextflow_schema.json`**
   - Update `peakcaller` parameter with new options

---

## 6. Configuration Parameters

### 6.1 New Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `condition_aware_norm` | Boolean | `false` | Enable condition-specific normalization |
| `normalisation_scope` | String | `global` | Scope for normalization: `global`, `condition`, `group`, `group_condition` |
| `pool_controls_by` | String | `condition` | How to pool IgG controls: `all`, `condition` |
| `gopeaks_pval` | Float | `0.05` | GoPeaks p-value threshold |
| `epic2_fdr` | Float | `0.01` | epic2 FDR threshold |
| `span_fdr` | Float | `0.05` | SPAN FDR threshold (default variant) |
| `span_stringent_fdr` | Float | `1.0E-9` | SPAN FDR threshold (stringent variant) |

### 6.2 Updated Parameters

| Parameter | Change |
|-----------|--------|
| `peakcaller` | Extended to support: `seacr`, `macs2_narrow`, `macs2_broad`, `gopeaks_narrow`, `gopeaks_broad`, `epic2_200bp`, `epic2_150bp`, `epic2_25bp`, `span_default`, `span_stringent` |
| `consensus_peak_mode` | Extended to support: `group`, `all`, `condition`, `group_condition` |

### 6.3 Parameter Schema Update

```json
{
  "peakcaller": {
    "type": "string",
    "default": "seacr",
    "description": "Comma-separated list of peak callers to run",
    "pattern": "^(seacr|macs2|macs2_narrow|macs2_broad|gopeaks|gopeaks_narrow|gopeaks_broad|epic2|epic2_200bp|epic2_150bp|epic2_25bp|span|span_default|span_stringent)(,(seacr|macs2|macs2_narrow|macs2_broad|gopeaks|gopeaks_narrow|gopeaks_broad|epic2|epic2_200bp|epic2_150bp|epic2_25bp|span|span_default|span_stringent))*$",
    "fa_icon": "fas fa-mountain"
  },
  "condition_aware_norm": {
    "type": "boolean",
    "default": false,
    "description": "Enable condition-specific normalization for differential analysis",
    "fa_icon": "fas fa-balance-scale"
  }
}
```

---

## 7. Output Structure

### 7.1 Extended Output Directory

```
results/
├── 01_prealign/
├── 02_alignment/
├── 03_peak_calling_preprocessing/
│   └── scale_factors/
│       ├── global_scale_factors.csv
│       └── condition_scale_factors.csv  # NEW
├── 04_called_peaks/
│   ├── seacr/
│   ├── macs2_narrow/                    # NEW (explicit)
│   ├── macs2_broad/                     # NEW
│   ├── gopeaks_narrow/                  # NEW
│   ├── gopeaks_broad/                   # NEW
│   ├── epic2_200bp/                     # NEW
│   ├── epic2_150bp/                     # NEW
│   ├── epic2_25bp/                      # NEW
│   ├── span_default/                    # NEW
│   ├── span_stringent/                  # NEW
│   ├── pooled_controls/                 # NEW
│   │   ├── WT_pooled_IgG.bam
│   │   └── KO_pooled_IgG.bam
│   └── consensus/
│       ├── by_group/
│       └── by_condition/                # NEW
├── 05_reporting/
│   ├── diffbind_samplesheets/           # NEW
│   │   ├── seacr_samplesheet.csv
│   │   ├── macs2_narrow_samplesheet.csv
│   │   ├── macs2_broad_samplesheet.csv
│   │   ├── gopeaks_narrow_samplesheet.csv
│   │   ├── epic2_200bp_samplesheet.csv
│   │   └── span_default_samplesheet.csv
│   └── peak_caller_comparison/          # NEW
│       ├── peak_counts_by_caller.tsv
│       ├── frip_by_caller.tsv
│       └── caller_overlap_heatmap.png
└── pipeline_info/
```

### 7.2 DiffBind Samplesheet Format

Generated samplesheets follow DiffBind requirements:

```csv
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
H3K36me2_WT_R1,CellLine,H3K36me2,WT,1,/path/to/H3K36me2_WT_R1.bam,IgG_WT,/path/to/WT_pooled_IgG.bam,/path/to/H3K36me2_WT_R1.narrowPeak,narrow
H3K36me2_WT_R2,CellLine,H3K36me2,WT,2,/path/to/H3K36me2_WT_R2.bam,IgG_WT,/path/to/WT_pooled_IgG.bam,/path/to/H3K36me2_WT_R2.narrowPeak,narrow
H3K36me2_KO_R1,CellLine,H3K36me2,KO,1,/path/to/H3K36me2_KO_R1.bam,IgG_KO,/path/to/KO_pooled_IgG.bam,/path/to/H3K36me2_KO_R1.narrowPeak,narrow
H3K36me2_KO_R2,CellLine,H3K36me2,KO,2,/path/to/H3K36me2_KO_R2.bam,IgG_KO,/path/to/KO_pooled_IgG.bam,/path/to/H3K36me2_KO_R2.narrowPeak,narrow
```

---

## 8. Testing Strategy

### 8.1 New Test Profiles

1. **`test_condition`**: Test condition-aware processing with minimal dataset
2. **`test_multi_caller`**: Test all peak callers on single sample
3. **`test_full_condition`**: Full experimental design with 2×2 (condition × replicate)

### 8.2 Test Samplesheet

**`conf/test_condition.config`:**
```groovy
params {
    config_profile_name        = 'Test condition-aware processing'
    config_profile_description = 'Test with condition metadata'

    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    input  = "${projectDir}/assets/samplesheet_test_condition.csv"
    genome = 'GRCh38'

    peakcaller = 'seacr,macs2_narrow,gopeaks_narrow'
    condition_aware_norm = true
    normalisation_scope = 'condition'
}
```

### 8.3 CI/CD Updates

Add to `.github/workflows/ci.yml`:
```yaml
- name: Run condition-aware test
  run: |
    nextflow run ${GITHUB_WORKSPACE} \
      -profile test_condition,docker \
      --outdir ./results
```

---

## 9. Migration Guide

### 9.1 Backward Compatibility

Existing samplesheets **without** the `condition` column will continue to work:
- `condition` defaults to `"default"` for all samples
- All processing remains global (no condition-specific handling)

### 9.2 Upgrading Samplesheets

To enable condition-aware processing:

1. Add `condition` column to samplesheet header
2. Populate condition values for all samples (including controls)
3. Set `--condition_aware_norm true` if needed
4. Optionally set `--normalisation_scope condition`

### 9.3 New Caller Usage

To use new peak callers:
```bash
# Single caller
nextflow run nf-core/cutandrun --peakcaller 'gopeaks_broad' ...

# Multiple callers
nextflow run nf-core/cutandrun --peakcaller 'seacr,macs2_narrow,epic2_200bp,span_default' ...

# All callers (for comparison)
nextflow run nf-core/cutandrun --peakcaller 'seacr,macs2_narrow,macs2_broad,gopeaks_narrow,gopeaks_broad,epic2_200bp,epic2_150bp,epic2_25bp,span_default,span_stringent' ...
```

---

## Appendix A: Resource Requirements

### A.1 Per-Caller Resources

| Caller | Memory | CPUs | Time (est.) | Notes |
|--------|--------|------|-------------|-------|
| SEACR | 4 GB | 2 | 5-15 min | Existing |
| MACS2 (narrow) | 4 GB | 1 | 5-15 min | |
| MACS2 (broad) | 4 GB | 1 | 5-15 min | |
| GoPeaks (narrow) | 4 GB | 1 | 5-10 min | |
| GoPeaks (broad) | 4 GB | 1 | 5-10 min | |
| epic2 (200bp) | 8 GB | 8 | 10-30 min | Multi-threaded |
| epic2 (150bp) | 8 GB | 8 | 10-30 min | Multi-threaded |
| epic2 (25bp) | 16 GB | 8 | 20-60 min | High memory |
| SPAN (default) | 8 GB | 1 | 15-45 min | JVM-based |
| SPAN (stringent) | 8 GB | 1 | 15-45 min | JVM-based |

### A.2 Base Config Labels

```groovy
// Add to conf/base.config
withLabel:process_epic2 {
    cpus   = { check_max( 8, 'cpus' ) }
    memory = { check_max( 8.GB * task.attempt, 'memory' ) }
    time   = { check_max( 2.h * task.attempt, 'time' ) }
}

withLabel:process_epic2_highres {
    cpus   = { check_max( 8, 'cpus' ) }
    memory = { check_max( 16.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h * task.attempt, 'time' ) }
}

withLabel:process_span {
    cpus   = { check_max( 1, 'cpus' ) }
    memory = { check_max( 8.GB * task.attempt, 'memory' ) }
    time   = { check_max( 2.h * task.attempt, 'time' ) }
}
```

---

## Appendix B: Genome Configuration

### B.1 Genome Size Mapping

| Species | Genome | MACS2 `--gsize` | epic2/SPAN `--genome` |
|---------|--------|-----------------|----------------------|
| Human | hg38/GRCh38 | `2.7E+9` | `hg38` |
| Human | hg19/GRCh37 | `2.7E+9` | `hg19` |
| Mouse | mm39/GRCm39 | `1.87E+9` | `mm10`* |
| Mouse | mm10/GRCm38 | `1.87E+9` | `mm10` |

*Note: mm39 auto-maps to mm10 for tools that don't recognize mm39.

### B.2 Chromosome Sizes Generation

For SPAN, generate chrom.sizes from BAM header:
```groovy
process GENERATE_CHROM_SIZES {
    input:
    tuple val(meta), path(bam)

    output:
    path("${meta.genome}.chrom.sizes"), emit: chrom_sizes

    script:
    """
    samtools view -H $bam | \\
        grep "^@SQ" | \\
        sed 's/@SQ\\tSN:\\([^\\t]*\\)\\tLN:\\([0-9]*\\).*/\\1\\t\\2/' > ${meta.genome}.chrom.sizes
    """
}
```

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2024-12-30 | Claude Code | Initial specification |
