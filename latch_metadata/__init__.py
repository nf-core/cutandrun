from latch.types.directory import LatchDir
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    NextflowMetadata,
    NextflowRuntimeResources,
    Params,
    Section,
    Spoiler,
    Text,
)

from .parameters import generated_parameters

flow = [
    Section(
        "Input",
        Params("input"),
        Spoiler(
            "Optional",
            Params("use_control"),
        ),
    ),
    Section(
        "Reference Genome Options",
        Fork(
            "genome_source",
            "",
            latch_genome_source=ForkBranch("Latch Certified Reference Genome", Params("latch_genome")),
            input_ref=ForkBranch(
                "Custom Reference Genome",
                Params("fasta", "gtf"),
                Spoiler(
                    "Reference Options",
                    Params("gene_bed", "save_reference"),
                    Text(
                        "The transcriptome and GTF files in iGenomes are vastly out of date with respect to current annotations from Ensembl e.g. human iGenomes annotations are from Ensembl release 75, while the current Ensembl release is 108. Please consider downloading and using a more updated version of your reference genome."
                    ),
                    Params("genome"),
                ),
            ),
        ),
    ),
    Spoiler(
        "Spikein Genome Options",
        Params(
            "spikein_genome",
            "spikein_bowtie2",
            "spikein_fasta",
            "save_spikein_aligned",
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Optional Parameters",
        Spoiler(
            "Alignment Filters",
            Params(
                "bowtie2",
                "blacklist",
                "remove_mitochondrial_reads",
                "mito_name",
                "minimum_alignment_q_score",
                "save_unaligned",
                "save_align_intermed",
                "remove_linear_duplicates",
                "end_to_end",
            ),
        ),
        Spoiler(
            "Peak Calling Options",
            Params("peakcaller"),
            Spoiler(
                "MACS2 Optional Parameters",
                Params(
                    "macs2_qvalue",
                    "macs2_pvalue",
                    "macs_gsize",
                    "macs2_narrow_peak",
                    "macs2_broad_cutoff",
                    "consensus_peak_mode",
                    "replicate_threshold",
                    "extend_fragments",
                ),
            ),
            Spoiler(
                "SEACR Options",
                Params(
                    "seacr_peak_threshold",
                    "seacr_norm",
                    "seacr_stringent",
                ),
            ),
        ),
        Spoiler(
            "Read trimming and QC Options",
            Params(
                "save_merged_fastq",
                "save_trimmed",
                "clip_r1",
                "clip_r2",
                "three_prime_clip_r1",
                "three_prime_clip_r2",
                "trim_nextseq",
                "dedup_target_reads",
            ),
        ),
        Spoiler(
            "Normalization Options",
            Params("normalisation_mode", "normalisation_binsize", "igg_scale_factor", "dump_scale_factors"),
        ),
        Spoiler(
            "Reporting/Plotting Options",
            Params(
                "dt_qc_bam_binsize",
                "dt_qc_corr_method",
                "dt_calc_all_matrix",
                "dt_heatmap_gene_beforelen",
                "dt_heatmap_gene_bodylen",
                "dt_heatmap_gene_afterlen",
                "dt_heatmap_peak_beforelen",
                "dt_heatmap_peak_afterlen",
                "min_frip_overlap",
                "min_peak_overlap",
                "igv_show_gene_names",
                "igv_sort_by_groups",
            ),
        ),
        Spoiler(
            "Flow Switching Options",
            Params(
                "only_input",
                "only_genome",
                "only_preqc",
                "only_alignment",
                "only_filtering",
                "only_peak_calling",
            ),
        ),
        Spoiler(
            "Skip Tools",
            Params(
                "skip_fastqc",
                "skip_trimming",
                "skip_removeduplicates",
                "skip_reporting",
                "skip_preseq",
                "skip_igv",
                "skip_dt_qc",
                "skip_heatmaps",
                "skip_peak_qc",
                "skip_multiqc",
            ),
        ),
        Spoiler(
            "MultiQC Options",
            Params(
                "multiqc_title",
                "multiqc_methods_description",
            ),
        ),
    ),
]


NextflowMetadata(
    display_name="nf-core/cutandrun",
    author=LatchAuthor(
        name="nf-core",
    ),
    repository="https://github.com/latchbio-nfcore/cutandrun",
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=100,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
    flow=flow,
)
