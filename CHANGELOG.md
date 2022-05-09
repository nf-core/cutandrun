# nf-core/cutandrun: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
8. Create bedGraph files ([`bedtools`](https://github.com/arq5x/bedtools2/)
9. Create bigWig coverage files ([`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
10. Peak calling specifically tailored for low background noise experiments ([`SEACR`](https://github.com/FredHutch/SEACR))
11. Consensus peak merging and reporting ([`bedtools`](https://github.com/arq5x/bedtools2/))
12. Quality control and analysis:
    1. Alignment, fragment length and peak analysis and replicate reproducibility ([`python`](https://www.python.org/))
    2. Heatmap peak analysis ([`deepTools`](https://github.com/deeptools/deepTools/))
13. Genome browser session ([`IGV`](https://software.broadinstitute.org/software/igv/))
14. Present QC for raw read, alignment and duplicate reads ([`MultiQC`](http://multiqc.info/))

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
- [[#75](https://github.com/nf-core/cutandrun/issues/75)] - Fixed error caused by emtpy peak files being passed to the `CALCULATE_FRIP` and `CALCULATE_PEAK_REPROD` python reporting modules
- [[#83]](https://github.com/nf-core/cutandrun/issues/83) - Fixed error in violin chart generation with cast to int64

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `samtools` | 1.13        | 1.14        |

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
| `picard`   | 2.25.7      | 2.27.1      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.
