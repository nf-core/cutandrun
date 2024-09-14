import typing
from dataclasses import dataclass
from enum import Enum

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchRule, Multiselect, NextflowParameter


@dataclass
class SampleSheet:
    group: str
    replicate: int
    fastq_1: LatchFile
    fastq_2: LatchFile
    control: str


class Reference(Enum):
    hg38 = "GRCh38 (Homo Sapiens hg38)"
    hg19 = "GRCh37 (Homo Sapiens hg19)"


generated_parameters = {
    "input": NextflowParameter(
        type=typing.List[SampleSheet],
        display_name="Input",
        description="Information about the samples in the experiment",
        section_title=None,
        samplesheet_type="csv",
        samplesheet=True,
        default=None,
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        default=None,
        display_name="Outdir",
        description="The output directory where the results will be saved. You have to use absolute paths to store on Cloud infrastructure.",
    ),
    "multiqc_title": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="MultiQC Title",
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "save_reference": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Reference",
        description="Save genome reference data to the output directory",
    ),
    "save_merged_fastq": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Merged FastQ",
        description="Save any technical replicate FASTQ files that were merged to the output directory",
    ),
    "genome_source": NextflowParameter(
        type=str,
        display_name="Reference Genome",
        description="Choose Reference Genome",
    ),
    "latch_genome": NextflowParameter(
        type=typing.Optional[Reference],
        default=None,
        section_title="Reference genome options",
        display_name="Genome",
        description="Name of the Latch Verfied Reference Genome",
    ),
    "save_trimmed": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Trimmed",
        description="Save trimmed FASTQ files to the output directory",
    ),
    "save_spikein_aligned": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Spikein Aligned",
        description="Save BAM files aligned to the spike-in genome to the output directory",
    ),
    "save_unaligned": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Unaligned",
        description="Save unaligned sequences to the output directory",
    ),
    "save_align_intermed": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Save Align Intermed",
        description="Save alignment intermediates to the output directory (WARNING: can be very large)",
    ),
    "email": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Email",
        description="Email address for completion summary.",
    ),
    "dump_scale_factors": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Dump Scale Factors",
        description="Output calculated scale factors from pipeline",
    ),
    "genome": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Genome",
        description="Name of iGenomes reference.",
    ),
    "bowtie2": NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        display_name="Bowtie2",
        description="Path to bowtie2 index",
    ),
    "gtf": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="GTF",
        description="Path to GTF annotation file",
    ),
    "gene_bed": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Gene Bed",
        description="Path to gene BED file",
    ),
    "blacklist": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Blacklist",
        description="Path to genome blacklist",
    ),
    "spikein_genome": NextflowParameter(
        type=typing.Optional[str],
        default="K12-MG1655",
        section_title=None,
        display_name="Spikein Genome",
        description="Name of the iGenome reference for the spike-in genome, defaulting to E. coli K12, for yeast set to R64-1-1, for fruit fly BDGP6",
    ),
    "spikein_bowtie2": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Spikein Bowtie2",
        description="Path to spike-in bowtie2 index",
    ),
    "spikein_fasta": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Spikein Fasta",
        description="Path to spike-in fasta",
    ),
    "fasta": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Fasta",
        description="Path to FASTA genome file.",
    ),
    "only_input": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only Input",
        description="Run pipeline up to input checking",
    ),
    "only_genome": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only Genome",
        description="Run pipeline up to reference preparation",
    ),
    "only_preqc": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only PreQC",
        description="Run pipeline up to pre-alignment",
    ),
    "only_alignment": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only Alignment",
        description="Run pipeline up to alignment",
    ),
    "only_filtering": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only Filtering",
        description="Run pipeline up to q-filtering",
    ),
    "only_peak_calling": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Only Peak Calling",
        description="Run pipeline up to peak calling",
    ),
    "skip_fastqc": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip FastQC",
        description="Skips fastqc reporting",
    ),
    "skip_trimming": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Trimming",
        description="Skips trimming",
    ),
    "skip_removeduplicates": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Remove Duplicates",
        description="Skips de-duplication",
    ),
    "skip_reporting": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Reporting",
        description="Skips reporting",
    ),
    "skip_preseq": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Preseq",
        description="Skips preseq reporting",
    ),
    "skip_igv": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip IGV",
        description="Skips igv session generation",
    ),
    "skip_dt_qc": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip DT QC",
        description="Skips deeptools QC reporting",
    ),
    "skip_heatmaps": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Heatmaps",
        description="Skips deeptools heatmap generation",
    ),
    "skip_peak_qc": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip Peak QC",
        description="Skips peak QC reporting",
    ),
    "skip_multiqc": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Skip MultiQC",
        description="Skips multiqc",
    ),
    "clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        display_name="Clip R1",
        description="Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).",
    ),
    "clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        display_name="Clip R2",
        description="Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        display_name="Three Prime Clip R1",
        description="Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        display_name="Three Prime Clip R2",
        description="Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.",
    ),
    "trim_nextseq": NextflowParameter(
        type=typing.Optional[int],
        default=0,
        section_title=None,
        display_name="Trim NextSeq",
        description="Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails.",
    ),
    "minimum_alignment_q_score": NextflowParameter(
        type=typing.Optional[int],
        default=20,
        section_title=None,
        display_name="Minimum Alignment Q Score",
        description="Filter reads below a q-score threshold",
    ),
    "remove_mitochondrial_reads": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Remove Mitochondrial Reads",
        description="Filter mitochondrial reads",
    ),
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "mito_name": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Mito Name",
        description="Name of mitochondrial reads in reference genome. Only necessary when using a custom (non-igenomes) reference genome.",
    ),
    "dedup_target_reads": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Dedup Target Reads",
        description="De-duplicate target reads AND control reads (default is control only)",
    ),
    "remove_linear_duplicates": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Remove Linear Duplicates",
        description="De-duplicate reads based on read 1 5' start position. Relevant for assays using linear amplification with tagmentation (default is false).",
    ),
    "end_to_end": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="End to End",
        description="Use --end-to-end mode of Bowtie2 during alignment",
    ),
    "normalisation_mode": NextflowParameter(
        type=typing.Optional[str],
        default="Spikein",
        section_title=None,
        display_name="Normalisation Mode",
        appearance_type=Multiselect(["Spikein", "RPKM", "CPM", "BPM", "None"], allow_custom=False),
        description='Sets the target read normalisation mode. Options are: ["Spikein", "RPKM", "CPM", "BPM", "None" ]',
    ),
    "normalisation_binsize": NextflowParameter(
        type=typing.Optional[int],
        default=50,
        section_title=None,
        display_name="Normalisation Binsize",
        description='If normalisation option is one of  "RPKM", "CPM", "BPM" - then the binsize that the reads count is calculated on is used.',
    ),
    "peakcaller": NextflowParameter(
        type=typing.Optional[str],
        default="seacr",
        section_title=None,
        display_name="Peakcaller",
        appearance_type=Multiselect(["seacr", "macs2"], allow_custom=False),
        description="Selects the peak caller for the pipeline. Options are: [seacr, macs2]. More than one peak caller can be chosen and the order specifies which is a primary peak called (the first) that will be used downstream. Any secondary peak callers will be run and outputed to the results folder.",
    ),
    "use_control": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="Use Control",
        description="Specifies whether to use a control to normalise peak calls against (e.g. IgG)",
    ),
    "extend_fragments": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="Extend Fragments",
        description="Specifies whether to extend paired-end fragments between the read mates when calculating coverage tracks",
    ),
    "igg_scale_factor": NextflowParameter(
        type=typing.Optional[float],
        default=0.5,
        section_title=None,
        display_name="IGG Scale Factor",
        description="Specifies whether the background control is scaled prior to being used to normalise peaks.",
    ),
    "seacr_peak_threshold": NextflowParameter(
        type=typing.Optional[float],
        default=0.05,
        section_title=None,
        display_name="SEACR Peak Threshold",
        description="SEACR specifies returns the top n fraction (between 0 and 1) of peaks based on total signal within peaks. This is only used if there are no controls included with the samples and if `--use_control` is `false`",
    ),
    "seacr_norm": NextflowParameter(
        type=typing.Optional[str],
        default="non",
        section_title=None,
        appearance_type=Multiselect(["non", "norm"], allow_custom=False),
        display_name="SEACR Norm",
        description="SEACR normalization.",
    ),
    "seacr_stringent": NextflowParameter(
        type=typing.Optional[str],
        default="stringent",
        section_title=None,
        appearance_type=Multiselect(["stringent", "relaxed"], allow_custom=False),
        display_name="SEACR Stringency",
        description="SEACR stringency.",
    ),
    "macs2_qvalue": NextflowParameter(
        type=typing.Optional[float],
        default=0.01,
        section_title=None,
        display_name="MACS2 Q-value",
        description="Q-value threshold for MACS2 peak caller.",
    ),
    "macs2_pvalue": NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        display_name="MACS2 P-value",
        description="P-value threshold for macs2 peak caller. If set it will override the q-value.",
    ),
    "macs_gsize": NextflowParameter(
        type=typing.Optional[float],
        default=2700000000.0,
        section_title=None,
        display_name="MACS Gsize",
        description="Parameter required by MACS2. If using an iGenomes reference these have been provided when `--genome` is set as *GRCh37*, *GRCh38*, *GRCm38*, *WBcel235*, *BDGP6*, *R64-1-1*, *EF2*, *hg38*, *hg19* and *mm10*. Otherwise the gsize will default to GRCh38.",
    ),
    "macs2_narrow_peak": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="MACS2 Narrow Peak",
        description="Determines whether MACS2 broad or narrow peak mode is used for the peak caller",
    ),
    "macs2_broad_cutoff": NextflowParameter(
        type=typing.Optional[float],
        default=0.1,
        section_title=None,
        display_name="MACS2 Broad Cutoff",
        description="MACS2 broad cutoff parameter",
    ),
    "consensus_peak_mode": NextflowParameter(
        type=typing.Optional[str],
        default="group",
        section_title=None,
        display_name="Consensus Peak Mode",
        appearance_type=Multiselect(["group", "all"], allow_custom=False),
        description="Specifies what samples to group together for consensus peaks. Options are [group, all]",
    ),
    "replicate_threshold": NextflowParameter(
        type=typing.Optional[float],
        default=1.0,
        section_title=None,
        display_name="Replicate Threshold",
        description="Minimum number of overlapping replicates needed for a consensus peak",
    ),
    "igv_show_gene_names": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="IGV Show Gene Names",
        description="Show gene names instead of symbols in IGV browser sessions",
    ),
    "dt_qc_bam_binsize": NextflowParameter(
        type=typing.Optional[int],
        default=500,
        section_title=None,
        display_name="DT QC BAM Binsize",
        description="Deeptools multiBamSummary bam bin size",
    ),
    "dt_qc_corr_method": NextflowParameter(
        type=typing.Optional[str],
        default="pearson",
        section_title=None,
        display_name="DT QC Corr Method",
        description="Deeptools Correlation Plot statistical calculation method",
    ),
    "dt_heatmap_gene_beforelen": NextflowParameter(
        type=typing.Optional[int],
        default=3000,
        section_title=None,
        display_name="DT Heatmap Gene Beforelen",
        description="Deeptools heatmap gene plot before length (bases)",
    ),
    "dt_heatmap_gene_bodylen": NextflowParameter(
        type=typing.Optional[int],
        default=5000,
        section_title=None,
        display_name="DT Heatmap Gene Bodylen",
        description="Deeptools heatmap gene plot body length (bases)",
    ),
    "dt_heatmap_gene_afterlen": NextflowParameter(
        type=typing.Optional[int],
        default=3000,
        section_title=None,
        display_name="DT Heatmap Gene Afterlen",
        description="Deeptools heatmap gene plot after length (bases)",
    ),
    "dt_heatmap_peak_beforelen": NextflowParameter(
        type=typing.Optional[int],
        default=3000,
        section_title=None,
        display_name="DT Heatmap Peak Beforelen",
        description="Deeptools heatmap peak plot before length (bases)",
    ),
    "dt_heatmap_peak_afterlen": NextflowParameter(
        type=typing.Optional[int],
        default=3000,
        section_title=None,
        display_name="DT Heatmap Peak Afterlen",
        description="Deeptools heatmap peak plot after length (bases)",
    ),
    "dt_calc_all_matrix": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="DT Calc All Matrix",
        description="Flag for whether to generate a heatmap for all samples together",
    ),
    "min_frip_overlap": NextflowParameter(
        type=typing.Optional[float],
        default=0.2,
        section_title=None,
        display_name="Min FRiP Overlap",
        description="Minimum fragment overlap for FriP score",
    ),
    "min_peak_overlap": NextflowParameter(
        type=typing.Optional[float],
        default=0.2,
        section_title=None,
        display_name="Min Peak Overlap",
        description="Minimum peak overlap for peak reproducibility plot",
    ),
    "igv_sort_by_groups": NextflowParameter(
        type=bool,
        default=True,
        section_title=None,
        display_name="IGV Sort by Groups",
        description="Sort the IGV output tracks by group",
    ),
    "singularity_pull_docker_container": NextflowParameter(
        type=bool,
        default=None,
        section_title=None,
        display_name="Singularity Pull Docker Container",
        description="Pull Docker container.",
    ),
    "multiqc_methods_description": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="MultiQC Methods Description",
        description="Custom MultiQC yaml file containing HTML including a methods description.",
    ),
}
