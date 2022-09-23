# ![nf-core/cutandrun](docs/images/nf-core-cutandrun_logo_light.png#gh-light-mode-only) ![nf-core/cutandrun](docs/images/nf-core-cutandrun_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/cutandrun/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/cutandrun/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/cutandrun/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/cutandrun/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/cutandrun/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.5653535-1073c8)](https://doi.org/10.5281/zenodo.5653535)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/cutandrun)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23cutandrun-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/cutandrun)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/cutandrun** is a best-practice bioinformatic analysis pipeline for CUT&RUN and CUT&Tag experimental protocols that were developed to study protein-DNA interactions and epigenomic profiling.

[CUT&RUN](https://elifesciences.org/articles/46314)

> Meers, M. P., Bryson, T. D., Henikoff, J. G., & Henikoff, S. (2019). Improved CUT&RUN chromatin profiling tools. _eLife_, _8_. https://doi.org/10.7554/eLife.46314

[CUT&Tag](https://www.nature.com/articles/s41467-019-09982-5)

> Kaya-Okur, H. S., Wu, S. J., Codomo, C. A., Pledger, E. S., Bryson, T. D., Henikoff, J. G., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. _Nature Communications_, _10_(1), 1930. https://doi.org/10.1038/s41467-019-09982-5]

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable, reproducible manner. It is capable of using containerisation and package management making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process, which makes it easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules).

The pipeline has been developed with continuous integration (CI) and test driven development (TDD) at its core. nf-core code and module linting as well as a battery of over 100 unit and integration tests run on pull request to the main repository and on release of the pipeline. On official release, automated CI tests run the pipeline on a full-sized dataset on AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/cutandrun/results).

![pipeline_diagram](docs/images/cutandrun-flow-diagram-v3.0.png)

## Pipeline summary

1. Check input files
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
5. Alignment to both target and spike-in genomes ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
6. Filter on quality, sort and index alignments ([`samtools`](https://sourceforge.net/projects/samtools/files/samtools/))
7. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/))
8. Create bedGraph files ([`bedtools`](https://github.com/arq5x/bedtools2/)
9. Create bigWig coverage files ([`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
10. Peak calling ([`SEACR`](https://github.com/FredHutch/SEACR), [`MACS2`](https://github.com/macs3-project/MACS))
11. Consensus peak merging and reporting ([`bedtools`](https://github.com/arq5x/bedtools2/))
12. Library complexity ([preseq]([Preseq | The Smith Lab](http://smithlabresearch.org/software/preseq)))
13. Fragment-based quality control ([`deepTools`](https://github.com/deeptools/deepTools/))
14. Peak-based quality control ([`bedtools`](https://github.com/arq5x/bedtools2/), custom python)
15. Heatmap peak analysis ([`deepTools`](https://github.com/deeptools/deepTools/))
16. Genome browser session ([`IGV`](https://software.broadinstitute.org/software/igv/))
17. Present all QC in web-based report ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/cutandrun -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   - Typical command for CUT&Run/CUT&Tag analysis:

   ```bash
   nextflow run nf-core/cutandrun --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

See [usage docs](https://nf-co.re/cutandrun/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/cutandrun pipeline comes with documentation about the pipeline [usage](https://nf-co.re/cutandrun/usage), [parameters](https://nf-co.re/cutandrun/parameters) and [output](https://nf-co.re/cutandrun/output).

## Credits

nf-core/cutandrun was originally written by Chris Cheshire ([@chris-cheshire](https://github.com/chris-cheshire)) and Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK.

The pipeline structure and parts of the downstream analysis were adapted from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html).

We thank Harshil Patel ([@drpatelh](https://github.com/drpatelh)) and everyone in the Luscombe Lab ([@luslab](https://github.com/luslab)) for their extensive assistance in the development of this pipeline.

![[The Francis Crick Institute](https://www.crick.ac.uk/)](docs/images/crick_logo.png)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#cutandrun` channel](https://nfcore.slack.com/channels/cutandrun) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/cutandrun for your analysis, please cite it using the following doi: [10.5281/zenodo.5653535](https://doi.org/10.5281/zenodo.5653535)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
