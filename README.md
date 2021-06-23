# ![nf-core/cutandrun](docs/images/nf-core-cutandrun_logo.png)

**Analysis pipeline for CUT&RUN and CUT&TAG experiments that includes sequencing QC, spike-in normalisation, IgG control normalisation, peak calling and downstream peak analysis.**.

[![GitHub Actions CI Status](https://github.com/nf-core/cutandrun/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/cutandrun/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/cutandrun/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/cutandrun/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/cutandrun.svg)](https://hub.docker.com/r/nfcore/cutandrun)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23cutandrun-4A154B?logo=slack)](https://nfcore.slack.com/channels/cutandrun)

## Introduction

**nf-core/cutandrun** is a bioinformatics best-practise analysis pipeline for CUT&Run and CUT&Tag sequencing data analysis to study protein-DNA interactions and epigenomic profiling.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/cutandrun/results).

## Pipeline summary

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Alignment to both target and spike-in genomes ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
5. Filter on quality, sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
6. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
7. Create bedGraph files ([`BEDTools`](https://github.com/arq5x/bedtools2/)
8. Create bigWig coverage files ([`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
9. Peak calling specifically tailored for low background noise ([`SEACR`](https://github.com/FredHutch/SEACR))
10. Quality control and analysis:
    1. Alignment, fragment length and peak analysis and replicate reproducibility ([`Python`](https://www.python.org/))
    2. Differential peak analysis ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    3. Heatmap peak analysis ([`deepTools`](https://github.com/deeptools/deepTools/))
11. Genome browser session ([`IGV`](https://software.broadinstitute.org/software/igv/))
12. Present QC for raw read, alignment and duplicate reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/cutandrun -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    * Typical command for CUT&Run/CUT&Tag analysis:

        ```bash
        nextflow run nf-core/cutandrun \
            -profile <docker/singularity/podman/conda/institute> \
            --input samplesheet.csv \
            --genome GRCh37
        ```

See [usage docs](https://nf-co.re/cutandrun/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/cutandrun pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/cutandrun/usage) and [output](https://nf-co.re/cutandrun/output).

## Credits

nf-core/cutandrun was originally written by Chris Cheshire ([@chris-cheshire](https://github.com/chris-cheshire)) and Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK.

The pipeline structure and parts of the downstream analysis were adapted from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html).

We thank Harshil Patel ([@drpatelh](https://github.com/drpatelh)) and everyone in the Luscombe Lab ([@luslab](https://github.com/luslab)) for their extensive assistance in the development of this pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#cutandrun` channel](https://nfcore.slack.com/channels/cutandrun) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/cutandrun for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
