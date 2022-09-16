# nf-core/cutandrun: Output

## Introduction

This document describes the output produced by the pipeline. We try to show both good and bad outputs for each section of the output reporting. This document can be used as a general guide for CUT&RUN analysis. Unless specified all outputs shown are from the MultiQC report generated at the end of the pipeline run.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [nf-core/cutandrun: Output](#nf-corecutandrun-output)
  - [Introduction](#introduction)
  - [Pipeline overview](#pipeline-overview)
  - [Preprocessing](#preprocessing)
    - [Samplesheet check](#samplesheet-check)
    - [Fastq merging](#fastq-merging)
    - [FastQC](#fastqc)
    - [TrimGalore](#trimgalore)
  - [Alignment](#alignment)
    - [Bowtie 2](#bowtie-2)
  - [Alignment post-processing](#alignment-post-processing)
    - [samtools](#samtools)
    - [picard MarkDuplicates/RemoveDuplicates](#picard-markduplicatesremoveduplicates)
  - [Peak Calling](#peak-calling)
    - [Bam to bedgraph](#bam-to-bedgraph)
    - [Clip bedfiles](#clip-bedfiles)
    - [Bed to bigwig](#bed-to-bigwig)
    - [SEACR peak calling](#seacr-peak-calling)
    - [BEDtools](#bedtools)
  - [Reporting](#reporting)
    - [Python reporting](#python-reporting)
    - [MultiQC](#multiqc)
    - [IGV](#igv)
    - [Deeptools](#deeptools)
  - [Workflow reporting and genomes](#workflow-reporting-and-genomes)
    - [Reference genome files](#reference-genome-files)
    - [Pipeline information](#pipeline-information)

## Preprocessing

### Samplesheet check

The first step of the pipeline is to verify the samplesheet structure and experimental design to ensure that it is valid.

### Fastq merging

If multiple libraries/runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage documentation](https://nf-co.re/rnaseq/usage#samplesheet-input) to see how to specify these samples in the input samplesheet.

<details markdown="1">
<summary>Output files</summary>

- `01_prealign/merged_fastq/`
  - `*.merged.fastq.gz`: If `--save_merged_fastq` is specified, concatenated FastQ files will be placed in this directory.

</details>

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `01_prealign/pretrim_fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after adapter and quality trimming please refer to the FastQC reports in the `01_prealign/trimgalore/fastqc/` directory.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

We perform FastQC on reads before and after trimming. The descriptions below apply to both reporting unless explicitly mentioned.

### Sequence Counts

This first FastQC report provides a first look at how many reads your FASTQ files contain. Predicted duplicates are also shown in black however, we recommend using the PICARD duplicate reports as they are more detailed and accurate.

As a general rule of thumb for an abundant epitope like H3K27Me3, we recommend no fewer than 5M aligned reads, the sequence counts must reflect this number plus any unaligned error reads. Assuming a max alignment error rate of 20%, we recommend a minimum read count of 6M raw reads. Less abundant epitopes may require deeper sequencing for a good signal.

Other things to look out for in this plot is the consistency of reads across your various groups and replicates. Large differences in the number of reads within biological replicates may be indicative of human error, problems with the wet-lab protocol or with the sequencing process.

![MultiQC - FastQC sequence counts plot](images/output/mqc_01_fastqc_sequence_counts.png)

#### Sequence Quality

FastqQC provides several reports that look at sequence quality from different views. The mean quality scores plot shows the sequence quality along the length of the read. It is normal to see some drop off in quality towards the end of the read, especially with long-read sequencing (150 bp). The plot should be green and with modern sequencing on good quality samples, should be > 30 throughout the majority of the read. There are many factors that affect the read quality but users can expect drops in average score as they move towards primary tissue or other tissue types that are difficult to extra good quality CUT&RUN data from.

![mqc_02_fastqc_per_base_sequence_quality](images/output/mqc_02_fastqc_per_base_sequence_quality.png)

The Per-sequence quality score report shows a different view of sequencing quality, showing the distribution of scores for each sample. This chart will peak where the majority of the reads are scored. In modern Illumina sequencing this should be curve at the end of the chart in an area > 30.

![f](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_03_fastqc_per_sequence_quality_scores.png)

The Per-base sequence quality report is not applicable to CUT&RUN data generally. Any discordant sequence content at the beginning of the reads are common phenomenon for CUT&RUN reads. Failing to pass the Per base sequence content does not mean your data failed. It can be due to Tn5 preferential binding or what you might be detecting is the 10-bp periodicity that shows up as a sawtooth pattern in the length distribution. If so, this is normal and will not affect alignment or peak calling.

The Per-sequence GC content report shows the distribution of the GC content of all reads in a sample. This should be centred around the average GC content % for your target organism.

![mqc_04_fastqc_per_sequence_gc_content](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_04_fastqc_per_sequence_gc_content.png)

An unusually shaped distribution such as one with dual peaks could indicate a contaminated library or some other kind of biased subset. In the image below, we see a signifiant batch effect between two groups of samples run on different days. The samples in red most likely are contaminated with DNA that has a different GC content to that of the target organism.

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_05_fastqc_per_sequence_gc_content.png)

#### Overrepresented sequences

FastQC provides three reports that focus on overrepresented sequences at the read level. All three must be looked at both before and after trimming to gain a detailed view of the types of sequence duplication that your samples contain.

A normal high-throughput library will contain a diverse set of sequences, with no individual sequence making up a tiny fraction of the whole. Finding that a single sequence is very overrepresented in the set either means that it is highly biologically significant, or indicates that the library is contaminated, or not as diverse as you expected.

The sequence duplication level plot shows percentages of the library that has a particular range of duplication counts. For example, the plot below shows that for some samples ~10% of the library has > 100 duplicated sequences. While not particularly informative on its own, if problems are identified in subsequent reports, it can help to identify the scale of the duplication problem.

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_06_fastqc_sequence_duplication_levels.png)

The second two reports must be read together both before and after trimming for maximum insight. The first report, overrepresented sequences shows the % of identified sequences in each library; the second, adapter content, shows a cumulative plot of adapter content at each base position. Due to the short insert length for CUT&RUN, 25 bp sequencing is possible for samples. When the sequencing length increases, more of the sequencing adapters will be sequenced. These sequences should be trimmed off therefore it is important to look at both of these plots after trimming to check that the adapter content has been removed and there are only small levels of overrepresented sequences. A clear adapter content plot after trimming but where there are still significant levels of overrepresented sequences could indicate something biologically significant or experimental error such as contamination.

**A 25 b.p CUT&Tag experiment with clear reports even before trimming**

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_07_fastqc_clear_overrep.png)

**A 150 b.p. CUT&RUN experiment with significant adapter content before trimming but that is clear after trimming**

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_08_fastqc_adapter_content.png)

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_10_fastqc_overrepresented_sequences.png)

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_07_fastqc_clear_overrep.png)

**CUT&RUN experiment that shows over represented sequences even after trimming all the adapter content away**

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_11_fastqc_overrepresented_sequences.png)

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_12_fastqc_adapter_content.png) 

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_13_fastqc_overrepresented_sequences.png)

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_09_fastqc_clear_adapter.png)

### TrimGalore

<details markdown="1">
<summary>Output files</summary>

- `01_prealign/trimgalore/`
  - `*.fq.gz`: If `--save_trimmed` is specified, FastQ files **after** adapter trimming will be placed in this directory.
  - `*_trimming_report.txt`: Log file generated by Trim Galore!.
- `01_prealign/trimgalore/fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for read 1 (_and read2 if paired-end_) **after** adapter trimming.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper tool around Cutadapt and FastQC to perform quality and adapter trimming on FastQ files. By default, Trim Galore! will automatically detect and trim the appropriate adapter sequence.

> **NB:** TrimGalore! will only run using multiple cores if you are able to use more than > 5 and > 6 CPUs for single- and paired-end data, respectively. The total cores available to TrimGalore! will also be capped at 4 (7 and 8 CPUs in total for single- and paired-end data, respectively) because there is no longer a run-time benefit. See [release notes](https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019) and [discussion whilst adding this logic to the nf-core/atacseq pipeline](https://github.com/nf-core/atacseq/pull/65).

**Typical plot showing many small sequences being trimmed from reads**

![MultiQC - cutadapt trimmed sequence length plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_14_fastqc_cutadapt_trimmed.png)

## Alignment

### Bowtie 2

<details markdown="1">
<summary>Output files</summary>

- `02_alignment/bowtie2`

</details>

Adapter-trimmed reads are mapped to the target and spike-in genomes using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). The pipeline will output the `.bam` files with index and samtools stats for only the final set by default. For example, the full pipeline will only output picard duplicates processed files as this is the final step before peak calling. If the pipeline is run with `--only_align`, then the `bam` files from the initial sorting and indexing will be copied to the output directory as the other steps are not run.

If `--save_align_intermed` is specified then all the `bam` files from all stages will be copied over to the output directory.

If `--save_spikein_aligned` is specified then the spike-in alignment files will also be published.

MultiQC shows several alignment-based reports however, the most important is the alignment score plot. A typical plot for a well-defined genome will look like the image below with high alignment scores and low levels of multi-mapped and unaligned sequences. Low levels of alignment can be due to a multitude of different factors but is generally a strong sign that the input library was of a poor quality. That being said, if the total number of aligned reads is still above the required level for the target epitope abundance then it doesnt not mean the sample has failed as there still may be enough information to answer the biological question asked.

![MultiQC - Bowtie2 paired-end mapping stats](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_15_bowtie2_pe.png)

The MultiQC report also includes a spike-in alignment report. This plot is important for deciding whether to normalise your samples using the spike-in genome or by another method.  The default mode in the pipeline is to normalise stacked reads before peak calling for epitope abundance using spike-in normalisation.

Traditionally, E. coli DNA is carried along with bacterially-produced enzymes that are used in CUT&RUN and CUT&Tag experiments and gets tagmented non-specifically during the reaction. The fraction of total reads that map to the E.coli genome depends on the yield of epitope-targeted CUT&Tag, and so depends on the number of cells used and the abundance of that epitope in chromatin. Since a constant amount of protein is added to the reactions and brings along a fixed amount of E. coli DNA, E. coli reads can be used to normalize epitope abundance in a set of experiments.

Since the introduction of these techniques there are several factors that have reduced the usefulness of this type of normalisation in certain experimental conditions. Firstly, many commercially available kits now have very low levels of E.coli DNA in them, which therefore requires users to spike-in their own DNA for normalisation which is not always done. Secondly the normalisation approach is dependant on the cell count between samples being constant, which in our experience is quite difficult to achieve especially in tissue samples.

The image below shows a typical plot for samples that have the correct amount of spike-in DNA for normalisation. The target samples have usually < 1% spike-in alignment (but still with > 1000 reads to reach above noise thresholds) while IgG should have the most as it has a large abundance (2-5%). In this example, the IgG will be brought inline with the other epitopes enabling proper peak calling using the IgG as a background. 

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_16_spikein_bowtie2_pe_plot.png)

If you see very low spike-in levels for all samples, it is likely your Tn5 had no residual e.coli DNA and that no additional spike-in DNA was added. In this case spike-in normalisation cannot be used and normalisation must be switched to read count no normalisation at all. 

If you see a strange distribution of spike-in DNA alignment that does not fit with your knowledge of the relative abundance of your IgG and target epitopes, this is indicative of a problem with the spike-in process and again, another normalisation option should be chosen.

In the plot below, it may initially look as though the spike-in distribution is too varied to be useful however, the larger IgG spike-in alignments correspond to the target samples with more sequencing depth and more spike-in alignments. The alignment counts are also include with epitope abundance therefore, these samples are actually good candidates for spike-in normalisation. 

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_17_bowtie2_pe_plot.png)

### Library Complexity

To estimate library complexity and identify potentially over-sequenced libraries (or libraries with a low information content) we run preseq. [Library Complexity](http://smithlabresearch.org/software/preseq/) estimates the complexity of a library, showing how many additional unique reads are sequenced for increasing total read count. A shallow curve indicates complexity saturation. The dashed line shows a perfectly complex library where total reads = unique reads.

The plot below shows a group of samples where the majority of unique molecules are accounted for by 50M. The total molecules detected stretches beyond 250M indicating the library is over sequenced and that an identical future experiment could be sequenced to a lower depth without loosing information.

![plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_18_preseq_plot.png)

## Alignment post-processing

### samtools

<details markdown="1">
<summary>Output files</summary>

- `aligner/bowtie2/intermediate/`
  - `.filtered.bam`: If `--publish_align_intermeds` is specified the original BAM file containing read alignments to the target genome will be placed in this directory.
  - `.filtered.bam.bai`: BAI file for BAM.
- `aligner/bowtie2/intermediate/samtools_stats`
  - `.filtered.bam.*stats`: various statistics regarding the BAM files.

</details>

BAM files are filtered for a minimum quality score using [SAMtools](http://samtools.sourceforge.net/).

### PICARD MarkDuplicates/RemoveDuplicates

<details markdown="1">
<summary>Output files</summary>

- `02_alignment/bowtie2/target/markdup/`
  - `.markdup.bam`: Coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM file and so will be saved by default in the results directory.
  - `.markdup.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory.
- `02_alignment/bowtie2/target/markdup/picard_metrics`
  - `.markdup.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

By default, the pipeline uses [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) to _mark_ the duplicate reads identified amongst the alignments to allow you to gauge the overall level of duplication in your samples.

If your data includes IgG controls, these will additionally be de-duplicated. It is not the normal protocol to de-duplicate the target reads, however, if this is required, use the `--dedup_target_reads true` switch.

The plot below shows a typical CUT&Tag experiment that has its PCR cycles optimised. We see a low level of non-optical duplication (from library amplification) in the target samples but more in the IgG samples as the reads in these samples derive from non-specific tagmentation in the CUT&Tag reactions.

![MultiQC - Picard MarkDuplicates metrics plot](/Users/cheshic/dev/repos/luslab/cutandrun/docs/images/output/mqc_19_picard_markduplicates.png)

High levels of duplication are not necessarily  a problem as long as they are consistent across biological replicates or other comparable groups. Given that the target samples are not de-duplicated by default, if the balance of duplicate reads is off when comparing two samples, it may lead to inaccurate peak calling and subsequent spurious signals. High levels of non-optical duplication are indicative of over-amplified samples.

## Fragment-based QC

## Peak Calling

### Bam to bedgraph

<details markdown="1">
<summary>Output files</summary>

- `03_peak_calling/01_bam_to_bedgraph`
  - `*.bedgraph`: bedgraph coverage file.

</details>

Converts bam files to the bedgraph format.

### Clip bedfiles

### Bed to bigwig

<details markdown="1">
<summary>Output files</summary>

- `03_peak_calling/03_bed_to_bigwig`
  - `*.bigWig`: bigWig coverage file.

</details>

The [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format is an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). This mitigates the need to load the much larger BAM files for data visualisation purposes which will be slower and result in memory issues. The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.

### SEACR peak calling

<details markdown="1">
<summary>Output files</summary>

- `03_peak_calling/04_called_peaks/`
  - `.peaks*.bed`: BED file containing peak coordinates and peak signal.

</details>

[SEACR](https://github.com/FredHutch/SEACR) is a peak caller for data with low background-noise, so is well suited to CUT&Run/CUT&Tag data. SEACR can take in IgG control bedGraph files in order to avoid calling peaks in regions of the experimental data for which the IgG control is enriched. If `--igg_control false` is specified, SEACR calls enriched regions in target data by selecting the top 5% of regions by AUC by default. This threshold can be overwritten using `--peak_threshold`.

![Python reporting - peaks reproduced](images/py_reproduced_peaks.png)

![Python reporting - aligned fragments within peaks](images/py_frags_in_peaks.png)

### Bedtools

<details markdown="1">
<summary>Output files</summary>

- `seacr/`
  - `{group}.consensus_peaks.pdf`: schematic showing which consensus peaks are shared across replicates within groups
  - `all_peaks.consensus_peaks.pdf`: schematic showing which consensus peaks are shared across all samples
- `seacr/consensus_peaks`
  - `{group}.consensus.peaks.bed`: BED containing consensus peaks for each group
  - `all_peaks.consensus.peaks.bed`: BED containing consensus peaks across all samples

</details>

The merge function from [BEDtools](https://github.com/arq5x/bedtools2) is used to merge replicate peaks of the same experimental group to create a consensus peak set. This can then optionally be filtered for consensus peaks contributed to be a threshold number of replicates using `--replicate_threshold`. Additionally, the same workflow is run merging across all samples.

![Peak calling - group consensus peak plot](images/consensus_peaks.png)
![Peak calling - group consensus peak plot](images/all_consensus_peaks.png)

## Reporting

### Python reporting

<details markdown="1">
<summary>Output files</summary>

- `04_reporting/qc/`
  - `report.pdf`: PDF report of all plots.
  - `*.png`: individual plots featured in the PDF report.
  - `*.csv`: corresponding data used to produce the plot.

</details>

Additional QC and analysis pertaining particularly to CUT&Run and CUT&Tag data are reported in this module. This report was adapted in python from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html).

![Python reporting - fragment length distribution](images/py_frag_hist.png)

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

<details markdown="1">
<summary>Output files</summary>

- `04_reporting/multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### IGV

<details markdown="1">
<summary>Output files</summary>

- `04_reporting/igv/`
  - `igv_session.xml`: IGV session.
  - `*.txt`: IGV input file configurations.

</details>

An IGV session file will be created at the end of the pipeline containing the normalised bigWig tracks, per-sample peaks, target genome fasta and annotation GTF. Once installed, open IGV, go to File > Open Session and select the igv_session.xml file for loading.

> **NB:** If you are not using an in-built genome provided by IGV you will need to load the annotation yourself e.g. in .gtf and/or .bed format.

### Deeptools

<details markdown="1">
<summary>Output files</summary>

- `04_reporting/heatmaps/<gene/peak>/`
  - `.plotHeatmap.pdf`: heatmap PDF.
  - `.computeMatrix.mat.gz`: heatmap matrix.
  - `*.mat.tab`: matrix and heatmap configs.

</details>

[deeptools](https://github.com/deeptools/deepTools/) sub-tools computeMatrix and plotHeatmap are used to assess the distribution of fragments around genes and peaks.

## Workflow reporting and genomes

### Reference genome files

<details markdown="1">
<summary>Output files</summary>

- `00_genome/target/`
  - `*.fa`: If the `--save_reference` parameter is provided then all of the genome reference files will be placed in this directory.
- `00_genome/target/index/`
  - `bowtie2`: Directory containing target Bowtie2 indices.
- `00_genome/spikein/index/`
  - `bowtie2`: Directory containing spike-in Bowtie2 indices.

</details>

A number of genome-specific files are generated by the pipeline because they are required for the downstream processing of the results. If the `--save_reference` parameter is provided then these will be saved in the `00_genome/` directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
