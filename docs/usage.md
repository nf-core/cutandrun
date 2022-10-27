# nf-core/cutandrun: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/cutandrun/usage](https://nf-co.re/cutandrun/usage)

**nf-core/cutandrun** is a best-practice bioinformatic analysis pipeline for CUT&RUN and CUT&Tag experimental protocols that where developed to study protein-DNA interactions and epigenomic profiling.

## Samplesheet input

You will need to create a samplesheet file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with the correct data structure as shown in the examples below.

```bash
--input <path to samplesheet file>
```

An example sample sheet structure is shown below. This defines two target experimental groups for the histone marks h3k27me3 and h3k4me3 with two biological replicates per group. Each antibody target also has an IgG control. The two IgG experiments are configured as biological replicates in the same group named `igg_ctrl`. They are assigned as controls to the two other groups using the last `control` column. If there are an equal number of replicates assigned to the samples from the control group as is the case below, the IgG controls will automatically be assigned to the same replicate number. If there is a mismatch then the first replicate of the control group will be assigned to all.

```bash
group,replicate,fastq_1,fastq_2,control
h3k27me3,1,READ1_FASTQ.gz,READ2_FASTQ.gz,igg_ctrl
h3k27me3,2,READ1_FASTQ.gz,READ2_FASTQ.gz,igg_ctrl
h3k4me3,1,READ1_FASTQ.gz,READ2_FASTQ.gz,igg_ctrl
h3k4me3,2,READ1_FASTQ.gz,READ2_FASTQ.gz,igg_ctrl
igg_ctrl,1,READ1_FASTQ.gz,READ2_FASTQ.gz,
igg_ctrl,2,READ1_FASTQ.gz,READ2_FASTQ.gz,
```

| Column      | Description                                                                                                 |
| ----------- | ----------------------------------------------------------------------------------------------------------- |
| `group`     | Group identifier for sample. This will be identical for replicate samples from the same experimental group. |
| `replicate` | Integer representing replicate number.                                                                      |
| `fastq_1`   | Full path to FastQ file for read 1. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `fastq_2`   | Full path to FastQ file for read 2. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz".   |
| `control`   | String representing the control group in the `group` column to which this replicate is assigned to.         |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/cutandrun --input samplesheet.csv  --outdir <OUTDIR> --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/cutandrun
```

## Pipeline Configuration Options

### Flow and Output Configuration

There are some options detailed on the parameters page that are prefixed with `save`, `skip` or `only`. These are flow control options that allow for saving additional output to the results directory, skipping unwanted portions of the pipeline or running the pipeline up to a certain point, which can be useful for testing.

### Genome Configuration

The easiest way to run the pipeline is by using one of the pre-configured genomes that reflect the available genomes at [iGenomes]([AWS iGenomes](https://ewels.github.io/AWS-iGenomes/)). Assign `genome` to one of the key words for iGenomes and all the available reference data will be automatically fetched. The pipeline uses the following reference data:

- Target genome FASTA

- Target genome Bowtie2 Index

- Target genome GTF

- Target genome BED (will be generated from the GTF if not supplied)

- Target genome Blacklist (blacklist files for major genomes are included in the pipeline)

- Spike-in genome FASTA

- Spike-in genome Bowtie2 Index

If the `genome` parameter is not supplied, the user must provide all the target genome data themselves (except the gene BED file). The default spike-in genome is e.coli given that this is the natural spike-in product of the protein production process. However, it is possible to spike-in different DNA during the experimental protocol and then set the `spikein_genome` to the target organism.

### Read Filtering and Duplication

After alignment using Bowtie2, mapped reads are filtered to remove those which do not pass a minimum quality threshold. This threshold can be changed using the `minimum_alignment_q_score` parameter.

CUT&RUN and CUT&Tag both integrate adapters into the vicinity of antibody-tethered enzymes, and the exact sites of integration are affected by the accessibility of surrounding DNA. Given these experimental parameters, it is expected that there are many fragments which share common starting and end positions; thus, such duplicates are generally valid but would be filtered out by de-duplication tools. However, there will be a fraction of fragments that are present due to PCR duplication that cannot be separated.

Control samples such as those from IgG datasets have relatively high duplication rates due to non-specific interactions with the genome; therefore, it is appropriate to remove duplicates from control samples.

The default for the pipeline therefore is to only run de-duplication on control samples. If it is suspected that there is a heavy fraction of PCR duplicates present in the primary samples then the parameter `dedup_target_reads` can be set using

`--dedup_target_reads`

### Read Normalisation

The default mode in the pipeline is to normalise stacked reads before peak calling for epitope abundance using spike-in normalisation.

Traditionally, E. coli DNA is carried along with bacterially-produced enzymes that are used in CUT&RUN and CUT&Tag experiments and gets tagmented non-specifically during the reaction. The fraction of total reads that map to the E.coli genome depends on the yield of epitope-targeted CUT&Tag, and so depends on the number of cells used and the abundance of that epitope in chromatin. Since a constant amount of protein is added to the reactions and brings along a fixed amount of E. coli DNA, E. coli reads can be used to normalize epitope abundance in a set of experiments.

Since the introduction of these techniques there are several factors that have reduced the usefulness of this type of normalisation in certain experimental conditions. Firstly, many commercially available kits now have very low levels of E.coli DNA in them, which therefore requires users to spike-in their own DNA for normalisation which is not always done. Secondly the normalisation approach is dependant on the cell count between samples being constant, which in our experience is quite difficult to achieve especially in tissue samples.

For these reasons we provide several other modes of normalisation based on read count; however, it should be noted that this form of normalisation is more simplistic and does not take into account epitope abundance. These normalisation modes are performed by [Deeptools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), some are more relevant than others to this type of data, we recommend using CPM with a bin size of 1 as a default.

| Mode      | Description                                                                                               |
| --------- | --------------------------------------------------------------------------------------------------------- |
| `Spikein` | The default mode which normalises by E. coli DNA.                                                         |
| `RPKM`    | Reads Per Kilobase per Million mapped reads. More relevant for transcript based assays.                   |
| `CPM`     | Counts Per Million mapped reads = number of reads per bin / number of mapped reads. Default bin size is 1 |
| `BPM`     | number of reads per bin / sum of all reads per bin (in millions),                                         |
| `None`    | Disables normalisation.                                                                                   |

Normalisation mode can be changed by the parameter `--normalisation_mode`.

### Peak Calling

This pipeline currently provides peak calling via `SEACR` or `MACS2` using the `peakcaller` parameter. If control samples are provided in the sample sheet by default they will be used to normalise the called peaks against non-specific background noise. Control normalisation can be disabled using `--use_control`. Additionally it may be necessary to scale control samples being used as background, especially when read count normalisation methods have been used at earlier stages in the pipeline. To scale the control samples before peak calling, change the `--igg_scale_factor` parameter to a number between 0-1. Multiple peak callers can be run by using comma separated values e.g. `--peakcaller SEACR,MACS2`, in this mode the primary peak caller is the first in the list and will be used for downstream processing; any additional peak callers will simply output to the results directory.

### Consensus Peaks

After peak calling, consensus peaks will be calculated based on merging peaks within the same groups. The number of replicates required for a valid peak can be changed using `replicate_threshold`. In some situations a user may which to call consensus peaks based on all samples, this can be configured by changing the `consensus_peak_mode` parameter from `group` to `all`.

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/cutandrun releases page](https://github.com/nf-core/cutandrun/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)

2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)

3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
