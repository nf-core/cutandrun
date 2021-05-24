# nf-core/cutandrun: Output

## Introduction

This document describes the output produced by the pipeline. Plots are taken from the python report which details summary details and analyses specific to CUT&Run/CUT&Tag data, and MultiQC report, which summarises results from tools used, at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Preprocessing](#preprocessing)
  * [ENA FTP](#ena-ftp) - Download FastQ files via SRA / ENA / GEO ids
  * [cat](#cat) - Merge re-sequenced FastQ files
  * [FastQC](#fastqc) - Raw read QC
  * [TrimGalore](#trimgalore) - Adapter and quality trimming
* [Alignment](#alignment)
  * [Bowtie 2](#bowtie-2) - Align reads to target and spike-in genomes
* [Alignment post-processing](#alignment-post-processing)
  * [SAMtools](#samtools) - Quality filter, sort and index alignments
  * [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking
* [Other steps](#other-steps)
  * [Calculate scale factor](#scale-factor) - Normalise between samples
  * [BEDTools and bedGraphToBigWig](#bedtools-and-bedgraphtobigwig) - Create bigWig coverage files
* [Peak calling](#peak-calling)
  * [SEACR](#seacr) - Peak calling for high signal-noise data
  * [Deeptools](#deeptools) - Analysis of peaks
* [Summary and quality control](#summary-and-quality-control)
  * [DESeq2](#deseq2) - PCA plot and differential peak analysis
  * [Python reporting](#python-reporting)
  * [MultiQC](#multiqc) - Present QC for raw reads, alignment, read counting and sample similarity  
  * [IGV](#igv) - Genome browser for viewing bigWigs in relation to genes
* [Workflow reporting and genomes](#workflow-reporting-and-genomes)
  * [Reference genome files](#reference-genome-files) - Saving reference genome indices/files
  * [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### ENA FTP

<details markdown="1">
<summary>Output files</summary>

* `public_data/`
  * `samplesheet.csv`: Auto-created samplesheet that can be used to run the pipeline.
  * `*.fastq.gz`: Paired-end/single-end reads downloaded from the ENA / SRA.
* `public_data/md5/`
  * `*.md5`: Files containing `md5` sum for FastQ files downloaded from the ENA / SRA.
* `public_data/runinfo/`
  * `*.runinfo.tsv`: Original metadata file downloaded from the ENA
  * `*.runinfo_ftp.tsv`: Re-formatted metadata file downloaded from the ENA

</details>

Please see the [usage documentation](https://nf-co.re/cutandrun/usage#direct-download-of-public-repository-data) for a list of supported public repository identifiers and how to provide them to the pipeline. The final sample information for all identifiers is obtained from the ENA which provides direct download links for FastQ files as well as their associated md5sums. If download links exist, the files will be downloaded in parallel by FTP otherwise they will NOT be downloaded. This is intentional because the tools such as `parallel-fastq-dump`, `fasterq-dump`, `prefetch` etc require pre-existing configuration files in the users home directory which makes automation tricky across different platforms and containerisation.

### cat

<details markdown="1">
<summary>Output files</summary>

* `fastq/`
  * `*.merged.fastq.gz`: If `--save_merged_fastq` is specified, concatenated FastQ files will be placed in this directory.

</details>

If multiple libraries/runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage documentation](https://nf-co.re/rnaseq/usage#samplesheet-input) to see how to specify these samples in the input samplesheet.

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after adapter and quality trimming please refer to the FastQC reports in the `trimgalore/fastqc/` directory.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

### TrimGalore

<details markdown="1">
<summary>Output files</summary>

* `trimgalore/`
  * `*.fq.gz`: If `--save_trimmed` is specified, FastQ files **after** adapter trimming will be placed in this directory.
  * `*_trimming_report.txt`: Log file generated by Trim Galore!.
* `trimgalore/fastqc/`
  * `*_fastqc.html`: FastQC report containing quality metrics for read 1 (*and read2 if paired-end*) **after** adapter trimming.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper tool around Cutadapt and FastQC to perform quality and adapter trimming on FastQ files. By default, Trim Galore! will automatically detect and trim the appropriate adapter sequence.

> **NB:** TrimGalore! will only run using multiple cores if you are able to use more than > 5 and > 6 CPUs for single- and paired-end data, respectively. The total cores available to TrimGalore! will also be capped at 4 (7 and 8 CPUs in total for single- and paired-end data, respectively) because there is no longer a run-time benefit. See [release notes](https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019) and [discussion whilst adding this logic to the nf-core/atacseq pipeline](https://github.com/nf-core/atacseq/pull/65).

![MultiQC - cutadapt trimmed sequence length plot](images/mqc_cutadapt_trimmed.png)

## Alignment

### Bowtie 2

<details markdown="1">
<summary>Output files</summary>

* `aligner/bowtie2/intermediate/`
  * `.bam`: If `--publish_align_intermeds` is specified the original BAM file containing read alignments to the target genome will be placed in this directory.
  * `.bam.bai`: BAI file for BAM.
* `aligner/bowtie2/intermediate/samtools_stats`
  * `.bam.*stats`: various statistics regarding the BAM files.
* `aligner/bowtie2/spikein/`
  * `.bam`: BAM file of reads aligned to the spike-in genome
  * `.bam.bai`: BAI file for spike-in BAM.
* `aligner/bowtie2/spikein/samtools_stats`
  * `.bam.*stats`: various statistics regarding the spike-in BAM files.

</details>

Adapter-trimmed reads are mapped to the target and spike-in genomes using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). A genome index is required to run Bowtie2 which is created automatically from the genome fasta input. By default, the only alignment files output are the quality filtered, marked and/or deduplicated alignment files. To output all alignment files including those directly from the aligner, set `--publish_align_intermed true`.

![MultiQC - Bowtie2 paired-end mapping stats](images/mqc_bowtie2_pe.png)

## Alignment post-processing

###  SAMtools

<details markdown="1">
<summary>Output files</summary>

* `aligner/bowtie2/intermediate/`
  * `.filtered.bam`: If `--publish_align_intermeds` is specified the original BAM file containing read alignments to the target genome will be placed in this directory.
  * `.filtered.bam.bai`: BAI file for BAM.
* `aligner/bowtie2/intermediate/samtools_stats`
  * `.filtered.bam.*stats`: various statistics regarding the BAM files.

</details>

BAM files are filtered for a minimum quality score of 0 using [SAMtools](http://samtools.sourceforge.net/). If duplicate marking and deduplication is skipped, then these will be the final alignment files under `aligner/bowtie2/`. Otherwise, these files can be included in the output by specifying `--publish_align_intermed true`.

### picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

* `aligner/bowtie2/`
  * `.markdup.bam`: Coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM file and so will be saved by default in the results directory.
  * `.markdup.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory.
* `aligner/bowtie2/picard_metrics`
  * `.markdup.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

By default, the pipeline uses [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) to *mark* the duplicate reads identified amongst the alignments to allow you to gauge the overall level of duplication in your samples. If your data includes IgG controls, these will additionally be deduplicated. You can skip this step via the `--skip_markduplicates` parameter. By default, this is the final processing step for the target BAM files and will appear in `aligner/bowtie2/`. However, if `--skip_markduplicates true` is set, this step will be skipped.

![MultiQC - Picard MarkDuplicates metrics plot](images/mqc_picard_markduplicates.png)

## Other steps

### BEDTools and bedGraphToBigWig

<details markdown="1">
<summary>Output files</summary>

* `ucsc/`
  * `*.bigWig`: bigWig coverage file.

</details>

The [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format is an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). This mitigates the need to load the much larger BAM files for data visualisation purposes which will be slower and result in memory issues. The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.

## Peak calling

### SEACR

<details markdown="1">
<summary>Output files</summary>

* `seacr/`
  * `.peaks*.bed`: BED file containing peak coordinates and peak signal.

</details>

[SEACR](https://github.com/FredHutch/SEACR) is a peak caller for data with low background-noise, so is well suited to CUT&Run/CUT&Tag data. SEACR can take in IgG control bedGraph files in order to avoid calling peaks in regions of the experimental data for which the IgG control is enriched. If `--igg_control false` is specified, SEACR calls enriched regions in target data by selecting the top 5% of regions by AUC by default. This threshold can be overwritten using `--peak_threshold`.

![Python reporting - peaks reproduced](images/py_reproduced_peaks.png)

![Python reporting - aligned fragments within peaks](images/py_frags_in_peaks.png)

### Deeptools

<details markdown="1">
<summary>Output files</summary>

* `deeptools/heatmaps/`
  * `.plotHeatmap.pdf`: heatmap PDF.
  * `.computeMatrix.mat.gz`: heatmap matrix.
  * `*.mat.tab`: matrix and heatmap configs.

</details>

[Deeptools](https://github.com/deeptools/deepTools/) sub-tools computeMatrix and plotHeatmap are used to assess the distribution of fragments around genes and peaks.

##  Summary and quality control

###  DESeq2

<details markdown="1">
<summary>Output files</summary>

* `deseq2_qc/`
  * `*.plots.pdf`: File containing PCA and hierarchical clustering plots.
  * `*.dds.RData`: File containing R `DESeqDataSet` object  generated
        by DESeq2, with either an rlog or vst `assay` storing the
        variance-stabilised data.
  * `*pca.vals.txt`: Matrix of values for the first 2 principal components.
  * `*sample.dists.txt`: Sample distance matrix.
  * `R_sessionInfo.log`: File containing information about R, the OS and attached or loaded packages.
* `deseq2_qc/size_factors/`
  * `*.txt`, `*.RData`: Files containing DESeq2 sizeFactors per sample.

</details>

[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) is one of the most commonly used software packages to perform differential expression analysis for RNA-seq datasets, but is used in this instance to compare between experimental CUT&Run datasets.

The script included in the pipeline uses DESeq2 to normalise read counts across all of the provided samples in order to create a PCA plot and a clustered heatmap showing pairwise Euclidean distances between the samples in the experiment. These help to show the similarity between groups of samples and can reveal batch effects and other potential issues with the experiment.

### Python reporting

<details markdown="1">
<summary>Output files</summary>

* `reports/`
  * `report.pdf`: PDF report of all plots.
  * `*.png`: individual plots featured in the PDF report.
  * `*.csv`: corresponding data used to produce the plot.

</details>

Additional QC and analysis pertaining particularly to CUT&Run and CUT&Tag data are reported in this module. This report was adapted in python from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html).

![Python reporting - fragment length distribution](images/py_frag_hist.png)

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

### IGV

<details markdown="1">
<summary>Output files</summary>

* `igv/`
  * `igv_session.xml`: IGV session.
  * `*.txt`: IGV input file configurations.

</details>

An IGV session file will be created at the end of the pipeline containing the normalised bigWig tracks, per-sample peaks, target genome fasta and annotation GTF. Once installed, open IGV, go to File > Open Session and select the igv_session.xml file for loading.

> **NB:** If you are not using an in-built genome provided by IGV you will need to load the annotation yourself e.g. in .gtf and/or .bed format.

## Workflow reporting and genomes

### Reference genome files

<details markdown="1">
<summary>Output files</summary>

* `genome/`
  * `*.fa`: If the `--save_reference` parameter is provided then all of the genome reference files will be placed in this directory.
* `genome/index/`
  * `bowtie2`: Directory containing target Bowtie2 indices.
* `genome/spikein_index/`
  * `bowtie2`: Directory containing spike-in Bowtie2 indices.

</details>

A number of genome-specific files are generated by the pipeline because they are required for the downstream processing of the results. If the `--save_reference` parameter is provided then these will be saved in the `genome/` directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
