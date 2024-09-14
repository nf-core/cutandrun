from typing import List, Optional

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import Reference, SampleSheet, initialize, nextflow_runtime


@workflow(metadata._nextflow_metadata)
def nf_nf_core_cutandrun(
    input: List[SampleSheet],
    multiqc_title: Optional[str],
    save_reference: bool,
    save_merged_fastq: bool,
    save_trimmed: bool,
    save_spikein_aligned: bool,
    save_unaligned: bool,
    save_align_intermed: bool,
    email: Optional[str],
    dump_scale_factors: bool,
    genome_source: str,
    genome: Optional[str],
    latch_genome: Optional[Reference],
    bowtie2: Optional[LatchDir],
    gtf: Optional[LatchFile],
    gene_bed: Optional[LatchFile],
    blacklist: Optional[LatchFile],
    spikein_bowtie2: Optional[LatchFile],
    spikein_fasta: Optional[LatchFile],
    fasta: Optional[LatchFile],
    only_input: bool,
    only_genome: bool,
    only_preqc: bool,
    only_alignment: bool,
    only_filtering: bool,
    only_peak_calling: bool,
    skip_fastqc: bool,
    skip_trimming: bool,
    skip_removeduplicates: bool,
    skip_reporting: bool,
    skip_preseq: bool,
    skip_igv: bool,
    skip_dt_qc: bool,
    skip_heatmaps: bool,
    skip_peak_qc: bool,
    skip_multiqc: bool,
    remove_mitochondrial_reads: bool,
    run_name: str,
    mito_name: Optional[str],
    dedup_target_reads: bool,
    remove_linear_duplicates: bool,
    macs2_pvalue: Optional[float],
    singularity_pull_docker_container: bool,
    multiqc_methods_description: Optional[str],
    spikein_genome: Optional[str] = "K12-MG1655",
    clip_r1: Optional[int] = 0,
    clip_r2: Optional[int] = 0,
    three_prime_clip_r1: Optional[int] = 0,
    three_prime_clip_r2: Optional[int] = 0,
    trim_nextseq: Optional[int] = 0,
    minimum_alignment_q_score: Optional[int] = 20,
    end_to_end: bool = True,
    normalisation_mode: Optional[str] = "Spikein",
    normalisation_binsize: Optional[int] = 50,
    peakcaller: Optional[str] = "seacr",
    use_control: bool = True,
    extend_fragments: bool = True,
    igg_scale_factor: Optional[float] = 0.5,
    seacr_peak_threshold: Optional[float] = 0.05,
    seacr_norm: Optional[str] = "non",
    seacr_stringent: Optional[str] = "stringent",
    macs2_qvalue: Optional[float] = 0.01,
    macs_gsize: Optional[float] = 2700000000.0,
    macs2_narrow_peak: bool = True,
    macs2_broad_cutoff: Optional[float] = 0.1,
    consensus_peak_mode: Optional[str] = "group",
    replicate_threshold: Optional[float] = 1.0,
    igv_show_gene_names: bool = True,
    dt_qc_bam_binsize: Optional[int] = 500,
    dt_qc_corr_method: Optional[str] = "pearson",
    dt_heatmap_gene_beforelen: Optional[int] = 3000,
    dt_heatmap_gene_bodylen: Optional[int] = 5000,
    dt_heatmap_gene_afterlen: Optional[int] = 3000,
    dt_heatmap_peak_beforelen: Optional[int] = 3000,
    dt_heatmap_peak_afterlen: Optional[int] = 3000,
    dt_calc_all_matrix: bool = True,
    min_frip_overlap: Optional[float] = 0.2,
    min_peak_overlap: Optional[float] = 0.2,
    igv_sort_by_groups: bool = True,
    outdir: LatchOutputDir = LatchOutputDir("latch:///CutandRun"),
) -> None:
    """
    nf-core/cutandrun is a best-practice bioinformatic analysis pipeline for CUT&RUN, CUT&Tag, and TIPseq experimental protocols that were developed to study protein-DNA interactions and epigenomic profiling.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.5653535-1073c8)](https://doi.org/10.5281/zenodo.5653535)

    **nf-core/cutandrun** is a best-practice bioinformatic analysis pipeline for CUT&RUN, CUT&Tag, and TIPseq experimental protocols that were developed to study protein-DNA interactions and epigenomic profiling.

    [CUT&RUN](https://elifesciences.org/articles/46314)

    > Meers, M. P., Bryson, T. D., Henikoff, J. G., & Henikoff, S. (2019). Improved CUT&RUN chromatin profiling tools. _eLife_, _8_. https://doi.org/10.7554/eLife.46314

    [CUT&Tag](https://www.nature.com/articles/s41467-019-09982-5)

    > Kaya-Okur, H. S., Wu, S. J., Codomo, C. A., Pledger, E. S., Bryson, T. D., Henikoff, J. G., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. _Nature Communications_, _10_(1), 1930. https://doi.org/10.1038/s41467-019-09982-5]

    [TIPseq](https://rupress.org/jcb/article/220/12/e202103078/212821)

    > Bartlett, D. A., Dileep, V., Handa, T., Ohkawa, Y., Kimura, H., Henikoff, S., & Gilbert, D. M. (2021). High-throughput single-cell epigenomic profiling by targeted insertion of promoters (TIP-seq). Journal of Cell Biology, 220(12), e202103078. https://doi.org/10.1083/jcb.202103078

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a portable, reproducible manner. It is capable of using containerisation and package management making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process, which makes it easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules).

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

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

    ## Usage

    First, prepare a samplesheet with your input data that looks as follows:

    `samplesheet.csv`:

    ```
    group,replicate,fastq_1,fastq_2,control
    h3k27me3,1,h3k27me3_rep1_r1.fastq.gz,h3k27me3_rep1_r2.fastq.gz,igg_ctrl
    h3k27me3,2,h3k27me3_rep2_r1.fastq.gz,h3k27me3_rep2_r2.fastq.gz,igg_ctrl
    igg_ctrl,1,igg_rep1_r1.fastq.gz,igg_rep1_r2.fastq.gz,
    igg_ctrl,2,igg_rep2_r1.fastq.gz,igg_rep2_r2.fastq.gz,
    ```

    Each row represents a pair of fastq files (paired end).

    ## Pipeline output

    To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/cutandrun/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/cutandrun/output).

    ## Credits

    nf-core/cutandrun was originally written by Chris Cheshire ([@chris-cheshire](https://github.com/chris-cheshire)) and Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK.

    The pipeline structure and parts of the downstream analysis were adapted from the original CUT&Tag analysis [protocol](https://yezhengstat.github.io/CUTTag_tutorial/) from the [Henikoff Lab](https://research.fredhutch.org/henikoff/en.html). The removal of duplicates arising from linear amplification (also known as T7 duplicates) in the TIPseq protocol was implemented as described in the original [TIPseq paper](https://rupress.org/jcb/article/220/12/e202103078/212821).

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

    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
        pvc_name=pvc_name,
        input=input,
        outdir=outdir,
        multiqc_title=multiqc_title,
        save_reference=save_reference,
        save_merged_fastq=save_merged_fastq,
        save_trimmed=save_trimmed,
        save_spikein_aligned=save_spikein_aligned,
        save_unaligned=save_unaligned,
        save_align_intermed=save_align_intermed,
        email=email,
        dump_scale_factors=dump_scale_factors,
        genome=genome,
        genome_source=genome_source,
        latch_genome=latch_genome,
        bowtie2=bowtie2,
        gtf=gtf,
        gene_bed=gene_bed,
        blacklist=blacklist,
        spikein_genome=spikein_genome,
        spikein_bowtie2=spikein_bowtie2,
        spikein_fasta=spikein_fasta,
        fasta=fasta,
        only_input=only_input,
        only_genome=only_genome,
        only_preqc=only_preqc,
        only_alignment=only_alignment,
        only_filtering=only_filtering,
        only_peak_calling=only_peak_calling,
        skip_fastqc=skip_fastqc,
        skip_trimming=skip_trimming,
        skip_removeduplicates=skip_removeduplicates,
        skip_reporting=skip_reporting,
        skip_preseq=skip_preseq,
        skip_igv=skip_igv,
        skip_dt_qc=skip_dt_qc,
        skip_heatmaps=skip_heatmaps,
        skip_peak_qc=skip_peak_qc,
        skip_multiqc=skip_multiqc,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        trim_nextseq=trim_nextseq,
        minimum_alignment_q_score=minimum_alignment_q_score,
        remove_mitochondrial_reads=remove_mitochondrial_reads,
        run_name=run_name,
        mito_name=mito_name,
        dedup_target_reads=dedup_target_reads,
        remove_linear_duplicates=remove_linear_duplicates,
        end_to_end=end_to_end,
        normalisation_mode=normalisation_mode,
        normalisation_binsize=normalisation_binsize,
        peakcaller=peakcaller,
        use_control=use_control,
        extend_fragments=extend_fragments,
        igg_scale_factor=igg_scale_factor,
        seacr_peak_threshold=seacr_peak_threshold,
        seacr_norm=seacr_norm,
        seacr_stringent=seacr_stringent,
        macs2_qvalue=macs2_qvalue,
        macs2_pvalue=macs2_pvalue,
        macs_gsize=macs_gsize,
        macs2_narrow_peak=macs2_narrow_peak,
        macs2_broad_cutoff=macs2_broad_cutoff,
        consensus_peak_mode=consensus_peak_mode,
        replicate_threshold=replicate_threshold,
        igv_show_gene_names=igv_show_gene_names,
        dt_qc_bam_binsize=dt_qc_bam_binsize,
        dt_qc_corr_method=dt_qc_corr_method,
        dt_heatmap_gene_beforelen=dt_heatmap_gene_beforelen,
        dt_heatmap_gene_bodylen=dt_heatmap_gene_bodylen,
        dt_heatmap_gene_afterlen=dt_heatmap_gene_afterlen,
        dt_heatmap_peak_beforelen=dt_heatmap_peak_beforelen,
        dt_heatmap_peak_afterlen=dt_heatmap_peak_afterlen,
        dt_calc_all_matrix=dt_calc_all_matrix,
        min_frip_overlap=min_frip_overlap,
        min_peak_overlap=min_peak_overlap,
        igv_sort_by_groups=igv_sort_by_groups,
        singularity_pull_docker_container=singularity_pull_docker_container,
        multiqc_methods_description=multiqc_methods_description,
    )


LaunchPlan(
    nf_nf_core_cutandrun,
    "Test Data",
    {
        "input": [
            SampleSheet(
                group="h3k27me3",
                replicate=1,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k27me3_rep1_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k27me3_rep1_r2.fastq.gz"),
                control="igg",
            ),
            SampleSheet(
                group="h3k27me3",
                replicate=2,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k27me3_rep2_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k27me3_rep2_r2.fastq.gz"),
                control="igg",
            ),
            SampleSheet(
                group="h3k4me3",
                replicate=1,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k4me3_rep1_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k4me3_rep1_r2.fastq.gz"),
                control="igg",
            ),
            SampleSheet(
                group="h3k4me3",
                replicate=2,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k4me3_rep2_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/h3k4me3_rep2_r2.fastq.gz"),
                control="igg",
            ),
            SampleSheet(
                group="igg",
                replicate=1,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/igg_rep1_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/igg_rep1_r2.fastq.gz"),
                control="",
            ),
            SampleSheet(
                group="igg",
                replicate=2,
                fastq_1=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/igg_rep2_r1.fastq.gz"),
                fastq_2=LatchFile("s3://latch-public/test-data/35729/cutandrun/test_data/igg_rep2_r2.fastq.gz"),
                control="",
            ),
        ],
        "blacklist": LatchFile("s3://latch-public/test-data/35729/cutandrun/hg38-blacklist.bed"),
        "spikein_bowtie2": LatchFile("s3://latch-public/test-data/35729/cutandrun/e_coli_U00096_3.tar.gz"),
        "spikein_fasta": LatchFile("s3://latch-public/test-data/35729/cutandrun/e_coli_U00096_3.fa"),
        "genome_source": "input_ref",
        "fasta": LatchFile("s3://latch-public/test-data/35729/cutandrun/hg38-chr20.fa"),
        "bowtie2": LatchDir("s3://latch-public/test-data/35729/cutandrun/hg38-chr20-bowtie2/"),
        "gtf": LatchFile("s3://latch-public/test-data/35729/cutandrun/hg38-chr20-genes.gtf"),
        "gene_bed": LatchFile("s3://latch-public/test-data/35729/cutandrun/hg38-chr20-genes.bed"),
        "run_name": "Test_Run",
        "peakcaller": "seacr",
        "mito_name": "chrM",
        "remove_linear_duplicates": True,
        "remove_mitochondrial_reads": True,
    },
)
