# nf-core/cutandrun: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.2.3] - 2025-08-22

### Fixes
- Fixed logic error in linear deduplication where the reads with the lowest MAPQ were kept instead.


## [3.2.2] - 2024-02-01

### Enhancements

- Updated pipeline template to nf-core/tools `2.12`.
- Added option to dump calculated scale factor values from the pipeline using `--dump_scale_factors`
- Added options to group IGV tracks by samplesheet group or by file type using `--igv_sort_by_groups`

### Fixes

- Fixed error that caused mismapping of IgG controls to their targets in certain samplesheet configurations.

## [3.2.1] - 2023-10-22

### Enhancements

- Updated pipeline template to nf-core/tools `2.10`.

### Fixes

- Fixed error that caused one random sample to run for markduplicates with certain versions of nextflow. Adding an explicit `collect()` on genome files has fixed this issue.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `fastqc`   | 11.9        | 12.1        |
| `picard`   | 3.0.0       | 3.1.0       |
| `multiqc`  | 1.15        | 1.17        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [3.2] - 2023-08-31

### Major Changes

- [[#189](https://github.com/nf-core/cutandrun/pull/189)] - Duplicates arising from linear amplification can be now removed by setting `--remove_linear_duplicates true`. `false` is default. [Linear amplification](https://doi.org/10.1186/1471-2164-4-19) is used in the [TIPseq protocol](https://doi.org/10.1083/jcb.202103078) in which genomic DNA is cut with Tn5 loaded with T7 promoter sequence that gets inserted in the cut DNA fragment. The T7 promoter sequence is then used to perform in vitro transcription to produce copies of the cut fragment. These duplicates are referred to as linear duplicates. Recent iterations of the CUT&Tag protocol, such as [nano-CUT&Tag](https://doi.org/10.1038/s41587-022-01535-4), have also been modified to include a linear amplification step. Credit to teemuronkko for this.
- [[#208](https://github.com/nf-core/cutandrun/issues/208)] - Updated the genome blacklists file to more accurate CUT&RUN specific regions rather than the old ChIP-Seq ENCODE blacklist. This should improve mapping rates and reduce spurious peaks. Credit to Adrija K for this. [[The CUT&RUN suspect list of problematic regions of the genome](https://doi.org/10.1186/s13059-023-03027-3)]

### Enhancements

- Updated pipeline template to nf-core/tools `2.8`.
- [[#189](https://github.com/nf-core/cutandrun/pull/189)] - Mitochondrial reads can be filtered before peak calling by setting `--remove_mitochondrial_reads true`. `false` is default. If using a custom reference genome, user can specify the string that is used to denote the mitochondrial reads in the reference using the `--mito_name` parameter.
- [[#189](https://github.com/nf-core/cutandrun/pull/189)] - The user can now specify explicitly if `end-to-end` vs `local` mode of Bowtie2 should be used by setting `--end_to_end` to `true` or `false`. `true` is default. In the `end-to-end` mode, all read characters are included when optiming an alignment. If the `local` mode is specified, Bowtie2 might exclude characters from one or both ends of the read to maximise alignment scores.
- [[#189](https://github.com/nf-core/cutandrun/pull/189)] - Added the name of the peak caller in the consensus peaks to make it clearer which peaks were used in the downstream reporting steps.
- [[#196](https://github.com/nf-core/cutandrun/pull/196)] - Extended documentation for most common alternative spike-in genomes, i.e. yeast and fruit fly. Credit to smoe for this.
- The Preseq module `lcextrap` was moved from `local` to `nf-core`
- Updated all nf-core modules to latest versions.

### Fixes

- Standardised channel structure for the nf-core Bowtie2 `align` module in the local `align_bowtie2` and `prepare_genome` subworkflows to prevent file errors.
- Fixed error caused by altered channel structure of the nf-core `bedtools_intersect` module.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `bedtools` | 1.13        | 1.14        |
| `multiqc`  | 2.30.0      | 2.31.0      |
| `samtools` | 1.16.1      | 1.17        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [3.1] - 2023-02-20

### Major Changes

- IgG controls will now be analysed by the deeptools QC subworkflow giving greater visibility on the quality of control samples.
- Updated the MACS2 default parameters to better process PA-Tn5/PA-Mnase based experiments. The new defaults use the q-value of `0.01` as the default cutoff in place of the p-value. The defaults have also been updated to keep duplicate reads in the peak finding process and also to shift the model to better account for nucleosome positioning `--nomodel --shift -75 --extsize 150 --keep-dup all`
- Deeptools plotHeatmap will now run for all samples as well as for singles. This can be disabled using the parameter `--dt_calc_all_matrix false`
- Bowtie2 default parameters have been updated to use the `--dovetail` option. After careful consideration and literature review, we have decided that overlapping mates can occur in CUT&RUN data and are still valid reads. This is also the agreed parameterisation in similar pipelines and also on the 4D nucleome portal.

### Enhancements

- Updated pipeline template to nf-core/tools `2.7.2`.
- Updated pipeline syntax to conform to new Nextflow version standards.
- Some locally defined subworkflows/modules have now been added to nf-core and re-imported as official modules/subworkflows.

### Fixes

- Fixed confusing config warnings that were being displayed on legitmate parameter configurations.
- Fixed deeptools correlation plots that were showing low levels of correlation even in test data by changing the plot to use Pearson correlation.
- Corrected the SEACR p-value parameter description.
- Fixed output of Picard mark/remove duplicate files so that the sorted, indexed bams for all files always output to the results folder.
- Spike-in genome processes and checks no longer run when the normalisation mode is set to something other than `SpikeIn`.
- Pipeline will now fail gracefully when single-end reads are detected.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.13        | 1.14        |
| `samtools` | 1.15.1      | 1.16.1      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [3.0] - 2022-09-26

### Major Changes

- Major rework of the pipeline internal flow structure. Metadata from processes (such as read counts) was previously annotated to a channel dictionary that was passed through the pipeline where various reporting processes could use the data. This was interacting with quite a few bugs in the Nextflow pipeline resume feature, causing lots of processes to rerun unnecessarily on resume. Any metadata generated in the pipeline is now written to files and passed where necessary to consuming reporting processes. This has drastically improved the number of processes that incorrectly rerun on resume.

- Re-organized the pipeline into clearer sections, breaking related processes into sub-workflows where possible. This is for better readability, but also to prepare the pipeline for the major upcoming nf-core feature of re-usable sub-workflows. as part of this rework, the pipeline now has distinct sections for fragment-based QC and peak-based QC.

- All reporting has been moved into MultiQC where possible. All PDF-based charting has been removed. Other PDF reports such as heatmaps and upset plots are still generated.

- We have listened to user comments that there is no guide on how to interpret the results from the pipeline. In response, we have revamped the documentation in the `output.md` document to describe the reporting in much more depth including good and bad examples of reporting output where possible.

- [[#140](https://github.com/nf-core/cutandrun/issues/140)] - IGV browser output has been reworked. We first fixed the performance issues with long load times by including the genome index into the session folder. IGV output now includes peaks from all peak callers used in pipeline, not just the primary one. Users can now select whether the gene track exported with the IGV session contains gene symbols or gene names. Several visual changes have been made to improve the default appearance and order of tracks.

- Added PreSeq library complexity reporting.

- Added full suite of fragment-based deepTools QC using the `multiBAMSummary` module. We generate three reporting from this fragment dataset: PCA, correlation and fingerprint plots. This has replaced our previous python implementation of sample correlation calculation.

- All coverage tracks generated from reads now extend reads to full fragment length by default. We feel this creates more realistic coverage tracks for CUT&RUN and improves the accuracy of other fragment-based reports.

### Enhancements

- Updated pipeline template to nf-core/tools `2.5.1`.
- [[#149](https://github.com/nf-core/cutandrun/issues/149)] - Pipeline will now use a blacklist file if provided to create an include list for the genome.
- The FRiP score is now calculated based on extended read fragments and not just mapped reads.
- [[#138](https://github.com/nf-core/cutandrun/issues/138)] - Better sample sheet error reporting.
- Gene bed files will now be automatically created from the GTF file if not supplied.
- The default minimum q-score for read quality has been changed from 0 to 20.
- [[#156](https://github.com/nf-core/cutandrun/issues/156)] SEACR has been better parameterized with dedicated config values for stringency and normalization. Credit to `CloXD` for this.
- deepTools heatmap generation has been better parameterized with dedicated config values for the gene and peak region settings.
- Consensus peak count reporting has been added to MultiQC.
- Reviewed and updated CI tests for better code coverage.
- Updated all nf-core modules to latest versions.

### Fixes

- Fixed some bugs in the passing of MACS2 peak data through the pipeline in v2.0. MACS2 peaks will now be correctly used and reporting on in the pipeline.
- [[#135](https://github.com/nf-core/cutandrun/issues/135)] - Removed many of the yellow warnings that were appearing in the pipeline to do with resource config options for processes that were not run.
- [[#137](https://github.com/nf-core/cutandrun/issues/137)] - Fixed the `workflow.OnComplete` error.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.12        | 1.13        |
| `picard`   | 2.27.2      | 2.27.4      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [2.0] - 2022-05-15

### Major Changes

- [[#53](https://github.com/nf-core/cutandrun/issues/53)] - Complete redesign of the samplesheet input system. Controls are no longer hard-coded to `igg` and the process of assigning controls to samples has been simplified.
- [[#93](https://github.com/nf-core/cutandrun/issues/93)], [[#73](https://github.com/nf-core/cutandrun/issues/73)] - Additional sample normalisation options have been added. In addition to normalising using detected spike-in DNA, there are now several options for normalising against read depth instead as well as skipping normalisation entirely.
- [[#62](https://github.com/nf-core/cutandrun/issues/62)] - Added MACS2 as an optional peak caller. Peak calling can now be altered using the `peakcaller` variable. Both peak callers can be run using `--peakcaller SEACR,MACS2`, the primary caller is the first item in the list and will be used downstream while the secondary will be run and outputted to the results folder.
- [[#101](https://github.com/nf-core/cutandrun/issues/101)] - `v1.1` ran consensus peak calling at both the group level and for all samples. This was causing performance issues for larger sample sets. There is now a new `consensus_peak_mode` parameter that defaults to `group`. Consensus will only be run on all samples if this is changed to `all`.

### Enhancements

- Updated pipeline template to nf-core/tools `2.3.2`.
- Upgraded pipeline to support the new nf-core module configuration system.
- More robust CI testing. Over 213 tests now before any code is merged with the main code base.
- More control over which parts of the pipeline run. Explicit skipping has been implemented for every section of the pipeline.
- Added options for scaling control data before it is used to call peaks. This is especially useful when using read depth normalisation as this can sometimes result in few peaks being called due to high background levels.
- Added support for Bowtie2 large indexes.
- IGV auto-session builder now supports `gff` and `fna` file extensions.
- Bowtie2 alignment has been altered to run in `--end-to-end` mode only if trimming is skipped. If trimming is activated then it will run in `--local` mode.
- [[#88](https://github.com/nf-core/cutandrun/issues/88)] - Many processes have been optimized for resource utilization. Users will especially notice that single thread processes will now only request 1 core rather than 2.
- [[#63](https://github.com/nf-core/cutandrun/issues/63)] - Custom containers for python reporting have now been condensed into a single container and added to BioConda.
- [[#76](https://github.com/nf-core/cutandrun/issues/76)] - Standardized python versions across reporting modules.

### Fixes

- [[#120](https://github.com/nf-core/cutandrun/issues/120)] - DeepTools compute matrix/heatmaps now only runs if there are peaks detected.
- [[#99](https://github.com/nf-core/cutandrun/issues/99)] - Large upset plots were causing process crashes. Upset plots will now fail gracefully if the number of samples in the consensus group is more than 10.
- [[#95](https://github.com/nf-core/cutandrun/issues/95)] - Fixed FRIP calculation performance issues and crashes.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `samtools` | 1.14        | 1.15.1      |
| `bowtie2`  | 2.4.2       | 2.4.4       |
| `picard`   | 2.25.7      | 2.27.2      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [1.1] - 2022-01-05

### Enhancements & fixes

- Updated pipeline template to nf-core/tools `2.2`
- [[#71](https://github.com/nf-core/cutandrun/issues/71)] - Bumped Nextflow version `21.04.0` -> `21.10.3`
- Added pipeline diagram to [[README](https://github.com/nf-core/cutandrun/blob/master/README.md)]
- Upgraded all modules (local and nf-core) to support the new versioning system
- The module `getchromsizes` was submitted to nf-core and moved from `local` to `nf-core`
- Added support for GFF files in IGV session generation
- [[#57](https://github.com/nf-core/cutandrun/issues/57), [#66](https://github.com/nf-core/cutandrun/issues/66)] - Upgraded version reporting in multiqc to support both software version by module and unique software versions. This improves detection of multi-version software usage in the pipeline
- [[#54](https://github.com/nf-core/cutandrun/issues/54)] - Fixed pipeline error where dots in sample ids inside the sample sheet were not correctly handled
- [[#75](https://github.com/nf-core/cutandrun/issues/75)] - Fixed error caused by empty peak files being passed to the `CALCULATE_FRIP` and `CALCULATE_PEAK_REPROD` python reporting modules
- [[#83]](https://github.com/nf-core/cutandrun/issues/83) - Fixed error in violin chart generation with cast to int64

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `samtools` | 1.13        | 1.14        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [1.0.0] - 2021-11-05

Initial release of nf-core/cutandrun, created with the [nf-core](https://nf-co.re/) template.

This pipeline is a best-practice bioinformatic analysis pipeline for CUT&Run and CUT&Tag experimental protocols that where developed to study protein-DNA interactions and epigenomic profiling.

nf-core/cutandrun was originally written by Chris Cheshire ([@chris-cheshire](https://github.com/chris-cheshire)) and Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK.

The pipeline structure and parts of the downstream analysis were adapted from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html).

We thank Harshil Patel ([@drpatelh](https://github.com/drpatelh)) and everyone in the Luscombe Lab ([@luslab](https://github.com/luslab)) for their extensive assistance in the development of this pipeline.

### Pipeline summary

1. Check input files
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
5. Alignment to both target and spike-in genomes ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
6. Filter on quality, sort and index alignments ([`samtools`](https://sourceforge.net/projects/samtools/files/samtools/))
7. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/))
8. Create bedGraph files ([`bedtools`](https://github.com/arq5x/bedtools2/))
9. Create bigWig coverage files ([`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
10. Peak calling specifically tailored for low background noise experiments ([`SEACR`](https://github.com/FredHutch/SEACR))
11. Consensus peak merging and reporting ([`bedtools`](https://github.com/arq5x/bedtools2/))
12. Quality control and analysis:
    1. Alignment, fragment length and peak analysis and replicate reproducibility ([`python`](https://www.python.org/))
    2. Heatmap peak analysis ([`deepTools`](https://github.com/deeptools/deepTools/))
13. Genome browser session ([`IGV`](https://software.broadinstitute.org/software/igv/))
14. Present QC for raw read, alignment and duplicate reads ([`MultiQC`](http://multiqc.info/))
