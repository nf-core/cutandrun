# nf-core/cutandrun: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/cutandrun/usage](https://nf-co.re/cutandrun/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple replicates

The `group` identifier is the same when you have multiple biological replicates from the same experimental group, just increment the `replicate` identifier appropriately. A special case for `group` is if you have non-specific IgG antibody control data that can be used for normalising your experimental CUT&Run (OR CUT&Tag) data. In this case, the `group` name for the IgG control data _must_ be set to `igg`. Below is an example for a single target group in triplicate, complemented by an IgG control duplicate:

```bash
group,replicate,fastq_1,fastq_2
target,1,H3K27me3_S1_L001_R1.fastq.gz,H3K27me3_S1_L001_R2.fastq.gz
target,2,H3K27me3_S2_L001_R1.fastq.gz,H3K27me3_S2_L001_R2.fastq.gz
target,3,H3K27me3_S3_L001_R1.fastq.gz,H3K27me3_S3_L001_R2.fastq.gz
igg,1,IGG_S1_L001_R1.fastq.gz,IGG_S1_L001_R2.fastq.gz
igg,2,IGG_S2_L001_R1.fastq.gz,IGG_S2_L001_R2.fastq.gz
```

There are 4 use-cases for various combinations of experimental and IgG control replicate numbers that are note-worthy:

* One-to-one:
Across all experimental groups and IgG control, the number of replicates are the same and numbered in a uniform fashion. In this case, IgG normalisation at the peak-calling level is done for matching replicate numbers.
* Many experimental replicates, one IgG replicate:
Each experimental replicate will be normalised by the given single IgG replicate.
* Equal numbers of replicates across experimental groups and more than one IgG replicate:
If the user puts forward groups with equal amounts of uniformly numbered replicates, but different to the multiple IgG replicates, then a warning will show informing the user that just the first IgG replicate will be used for normalisation.
* Varying experimental and IgG replicates:
If the user puts forward groups with varying, ununiform experimental replicate numbers and more than one IgG replicate, an error will be thrown up and the pipeline will hault. If you wish to merge the IgG control data in this case, see next steps.

It is _recommended_ to have an IgG control for normalising your experimental data and this is the default action for the pipeline. However, if you run the pipeline without IgG control data you must supply `--igg_control false`.

### Multiple runs of the same library

The `group` and `replicate` identifiers are the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth), or if you would like to merge technical replicates. The pipeline will concatenate the raw reads before alignment. Below is an example for two samples, one experimental and one control, sequenced across multiple lanes:

```bash
group,replicate,fastq_1,fastq_2
target,1,H3K27me3_S1_L001_R1.fastq.gz,H3K27me3_S1_L001_R2.fastq.gz
target,1,H3K27me3_S1_L002_R1.fastq.gz,H3K27me3_S1_L002_R2.fastq.gz
igg,1,IGG_S1_L001_R1.fastq.gz,IGG_S1_L001_R2.fastq.gz
igg,1,IGG_S1_L002_R1.fastq.gz,IGG_S1_L002_R2.fastq.gz
```

### Full design

A final design file may look something like the one below. This is for one experimental group in triplicate where the last replicate of the `treatment` group has been sequenced twice, another experimental group in duplicate, and one IgG control group.

```bash
group,replicate,fastq_1,fastq_2
h3k27me3,1,H3K27me3_S1_L001_R1.fastq.gz,H3K27me3_S1_L001_R2.fastq.gz
h3k27me3,2,H3K27me3_S2_L001_R1.fastq.gz,H3K27me3_S2_L001_R2.fastq.gz
h3k27me3,3,H3K27me3_S3_L001_R1.fastq.gz,H3K27me3_S3_L001_R2.fastq.gz
h3k27me3,3,H3K27me3_S3_L002_R1.fastq.gz,H3K27me3_S3_L002_R2.fastq.gz
h3k4me3,1,H3K4me3_S1_L001_R1.fastq.gz,H3K4me3_S1_L001_R2.fastq.gz
h3k4me3,2,H3K4me3_S2_L001_R1.fastq.gz,H3K4me3_S2_L001_R2.fastq.gz
igg,1,IGG_S1_L001_R1.fastq.gz,IGG_S1_L001_R2.fastq.gz
igg,2,IGG_S2_L001_R1.fastq.gz,IGG_S2_L001_R2.fastq.gz
```

| Column         | Description                                                                                                 |
|----------------|-------------------------------------------------------------------------------------------------------------|
| `group`        | Group identifier for sample. This will be identical for replicate samples from the same experimental group. |
| `replicate`    | Integer representing replicate number.                                                                      |
| `fastq_1`      | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `fastq_2`      | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Direct download of public repository data

> **NB:** This is an experimental feature but should work beautifully when it does! :)

The pipeline has been set-up to automatically download and process the raw FastQ files from public repositories. Identifiers can be provided in a file, one-per-line via the `--public_data_ids` parameter. Currently, the following identifiers are supported:

| `SRA`        | `ENA`        | `GEO`      |
|--------------|--------------|------------|
| SRR11605097  | ERR4007730   | GSM4432381 |
| SRX8171613   | ERX4009132   | GSE147507  |
| SRS6531847   | ERS4399630   |            |
| SAMN14689442 | SAMEA6638373 |            |
| SRP256957    | ERP120836    |            |
| SRA1068758   | ERA2420837   |            |
| PRJNA625551  | PRJEB37513   |            |

If `SRR`/`ERR` run ids are provided then these will be resolved back to their appropriate `SRX`/`ERX` ids to be able to merge multiple runs from the same experiment. This is conceptually the same as merging multiple libraries sequenced from the same sample.

The final sample information for all identifiers is obtained from the ENA which provides direct download links for FastQ files as well as their associated md5 sums. If download links exist, the files will be downloaded in parallel by FTP otherwise they will NOT be downloaded. This is intentional because the tools such as `parallel-fastq-dump`, `fasterq-dump`, `prefetch` etc require pre-existing configuration files in the users home directory which makes automation tricky across different platforms and containerisation.

As a bonus, the pipeline will also generate a valid samplesheet with paths to the downloaded data that can be used with the `--input` parameter, however, it is highly recommended that you double-check that all of the identifiers you defined using `--public_data_ids` are represented in the samplesheet. Also, public databases don't reliably hold information such as strandedness information so you may need to amend these entries too. All of the sample metadata obtained from the ENA has been appended as additional columns to help you manually curate the samplesheet before you run the pipeline.

If you have a GEO accession (found in the data availability section of published papers) you can directly download a text file containing the appropriate SRA ids to pass to the pipeline:

* Search for your GEO accession on [GEO](https://www.ncbi.nlm.nih.gov/geo)
* Click `SRA Run Selector` at the bottom of the GEO accession page
* Select the desired samples in the `SRA Run Selector` and then download the `Accession List`

This downloads a text file called `SRR_Acc_List.txt` which can be directly provided to the pipeline e.g. `--public_data_ids SRR_Acc_List.txt`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/cutandrun --input samplesheet.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/cutandrun
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/cutandrun releases page](https://github.com/nf-core/cutandrun/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
    * Pulls software from Docker Hub: [`nfcore/cutandrun`](https://hub.docker.com/r/nfcore/cutandrun/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
    * Pulls software from Docker Hub: [`nfcore/cutandrun`](https://hub.docker.com/r/nfcore/cutandrun/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
    * Pulls software from Docker Hub: [`nfcore/cutandrun`](https://hub.docker.com/r/nfcore/cutandrun/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
    * Pulls software from Docker Hub: [`nfcore/cutandrun`](https://hub.docker.com/r/nfcore/cutandrun/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
    * Pulls software from Docker Hub: [`nfcore/cutandrun`](https://hub.docker.com/r/nfcore/cutandrun/)
* `conda`
    * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
    withName: star {
        memory = 32.GB
    }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
