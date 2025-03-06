/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

// Validate input parameters in specialised library
WorkflowCutandrun.initialise(params, log)
def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.blacklist,
    params.bowtie2,
    params.fasta,
    params.gtf,
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if(params.normalisation_mode == "Spikein") {
    // Check spike-in only if it is enabled
    checkPathParamList = [
        params.spikein_bowtie2,
        params.spikein_fasta
    ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
}

// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

ch_blacklist = Channel.empty()
if (params.blacklist) {
    ch_blacklist = Channel.from( file(params.blacklist) )
}
else {
    ch_blacklist = Channel.empty()
    WorkflowCutandrun.blacklistWarn(log)
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Stage awk files for parsing log files
ch_bt2_to_csv_awk     = file("$projectDir/bin/bt2_report_to_csv.awk"    , checkIfExists: true)
ch_dt_frag_to_csv_awk = file("$projectDir/bin/dt_frag_report_to_csv.awk", checkIfExists: true)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

// Load up and check multiqc base config and custom configs
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_frag_len_header_multiqc              = file("$projectDir/assets/multiqc/frag_len_header.txt", checkIfExists: true)
ch_frip_score_header_multiqc            = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_counts_header_multiqc           = file("$projectDir/assets/multiqc/peak_counts_header.txt", checkIfExists: true)
ch_peak_counts_consensus_header_multiqc = file("$projectDir/assets/multiqc/peak_counts_consensus_header.txt", checkIfExists: true)
ch_peak_reprod_header_multiqc           = file("$projectDir/assets/multiqc/peak_reprod_header.txt", checkIfExists: true)
ch_linear_duplication_header_multiqc    = file("$projectDir/assets/multiqc/linear_duplication_header.txt", checkIfExists: true)


/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

// Init aligners
def prepare_tool_indices = ["bowtie2"]

// Check peak caller params
def caller_list = ['seacr', 'macs2']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : ['seacr']
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { INPUT_CHECK                } from "../subworkflows/local/input_check"
include { CUT as PEAK_TO_BED         } from '../modules/local/linux/cut'
include { AWK as AWK_NAME_PEAK_BED   } from "../modules/local/linux/awk"
include { IGV_SESSION                } from "../modules/local/python/igv_session"
include { AWK as AWK_EXTRACT_SUMMITS } from "../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW        } from "../modules/local/samtools_custom_view"
include { FRAG_LEN_HIST              } from "../modules/local/python/frag_len_hist"
include { MULTIQC                    } from "../modules/local/multiqc"

/*
 * SUBWORKFLOWS
 */
include { PREPARE_GENOME                                   } from "../subworkflows/local/prepare_genome"
include { FASTQC_TRIMGALORE                                } from "../subworkflows/local/fastqc_trimgalore"
include { ALIGN_BOWTIE2                                    } from "../subworkflows/local/align_bowtie2"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_TARGET_META  } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_SPIKEIN_META } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_PICARD_DUP_META  } from "../subworkflows/local/extract_metadata_awk"
include { MARK_DUPLICATES_PICARD                           } from "../subworkflows/local/mark_duplicates_picard"
include { MARK_DUPLICATES_PICARD as DEDUPLICATE_PICARD     } from "../subworkflows/local/mark_duplicates_picard"
include { CONSENSUS_PEAKS                                  } from "../subworkflows/local/consensus_peaks"
include { CONSENSUS_PEAKS as CONSENSUS_PEAKS_ALL           } from "../subworkflows/local/consensus_peaks"
include { EXTRACT_FRAGMENTS                                } from "../subworkflows/local/extract_fragments"
include { PREPARE_PEAKCALLING                              } from "../subworkflows/local/prepare_peakcalling"
include { DEEPTOOLS_QC                                     } from "../subworkflows/local/deeptools_qc"
include { PEAK_QC                                          } from "../subworkflows/local/peak_qc"
include { SAMTOOLS_VIEW_SORT_STATS as FILTER_READS         } from "../subworkflows/local/samtools_view_sort_stats"
include { DEDUPLICATE_LINEAR                               } from "../subworkflows/local/deduplicate_linear"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CAT_FASTQ                                                    } from "../modules/nf-core/cat/fastq/main"
include { PRESEQ_LCEXTRAP                                              } from "../modules/nf-core/preseq/lcextrap/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                         } from "../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                       } from "../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG                         } from "../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG                       } from "../modules/nf-core/macs2/callpeak/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE      } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS     } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE          } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS         } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE_ALL  } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS_ALL } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE_ALL      } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS_ALL     } from "../modules/nf-core/deeptools/plotheatmap/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                                  } from "../modules/local/custom_dumpsoftwareversions"

/*
 * SUBWORKFLOWS
 */


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CUTANDRUN {

    // Init
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    if(params.run_genome_prep) {
        PREPARE_GENOME (
            prepare_tool_indices,
            ch_blacklist
        )
        ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.versions)
    }

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    if(params.run_input_check) {
        INPUT_CHECK (
            ch_input
        )

        INPUT_CHECK.out.reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split("_")[0..-2].join("_")
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
    }

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    if(params.run_cat_fastq) {
        CAT_FASTQ (
            ch_fastq.multiple
        )
        ch_software_versions = ch_software_versions.mix(CAT_FASTQ.out.versions)

        CAT_FASTQ.out.reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //ch_cat_fastq | view

    /*
     * SUBWORKFLOW: Read QC, trim adapters and perform post-trim read QC
     */
    if(params.run_trim_galore_fastqc) {
        FASTQC_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc,
            params.skip_trimming
        )
        ch_trimmed_reads     = FASTQC_TRIMGALORE.out.reads
        ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //FASTQC_TRIMGALORE.out.reads | view

    /*
    * SUBWORKFLOW: Alignment to target and spikein genome using botwtie2
    */
    ch_orig_bam                   = Channel.empty()
    ch_orig_spikein_bam           = Channel.empty()
    ch_bowtie2_log                = Channel.empty()
    ch_bowtie2_spikein_log        = Channel.empty()
    ch_samtools_bam               = Channel.empty()
    ch_samtools_bai               = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_samtools_spikein_bam       = Channel.empty()
    ch_samtools_spikein_bai       = Channel.empty()
    ch_samtools_spikein_stats     = Channel.empty()
    ch_samtools_spikein_flagstat  = Channel.empty()
    ch_samtools_spikein_idxstats  = Channel.empty()
    if(params.run_alignment) {
        if (params.aligner == "bowtie2") {
            ALIGN_BOWTIE2 (
                ch_trimmed_reads,
                PREPARE_GENOME.out.bowtie2_index,
                PREPARE_GENOME.out.bowtie2_spikein_index,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.spikein_fasta
            )
            ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.versions)
            ch_orig_bam                   = ALIGN_BOWTIE2.out.orig_bam
            ch_orig_spikein_bam           = ALIGN_BOWTIE2.out.orig_spikein_bam
            ch_bowtie2_log                = ALIGN_BOWTIE2.out.bowtie2_log
            ch_bowtie2_spikein_log        = ALIGN_BOWTIE2.out.bowtie2_spikein_log

            ch_samtools_bam               = ALIGN_BOWTIE2.out.bam
            ch_samtools_bai               = ALIGN_BOWTIE2.out.bai
            ch_samtools_stats             = ALIGN_BOWTIE2.out.stats
            ch_samtools_flagstat          = ALIGN_BOWTIE2.out.flagstat
            ch_samtools_idxstats          = ALIGN_BOWTIE2.out.idxstats

            ch_samtools_spikein_bam       = ALIGN_BOWTIE2.out.spikein_bam
            ch_samtools_spikein_bai       = ALIGN_BOWTIE2.out.spikein_bai
            ch_samtools_spikein_stats     = ALIGN_BOWTIE2.out.spikein_stats
            ch_samtools_spikein_flagstat  = ALIGN_BOWTIE2.out.spikein_flagstat
            ch_samtools_spikein_idxstats  = ALIGN_BOWTIE2.out.spikein_idxstats
        }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: extract aligner metadata
     */
    ch_metadata_bt2_target  = Channel.empty()
    ch_metadata_bt2_spikein = Channel.empty()
    if (params.aligner == "bowtie2" && params.run_alignment) {
        EXTRACT_BT2_TARGET_META (
            ch_bowtie2_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_target = EXTRACT_BT2_TARGET_META.out.metadata
        ch_software_versions   = ch_software_versions.mix(EXTRACT_BT2_TARGET_META.out.versions)

        EXTRACT_BT2_SPIKEIN_META (
            ch_bowtie2_spikein_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_spikein = EXTRACT_BT2_SPIKEIN_META.out.metadata
    }
    //ch_metadata_bt2_target | view
    //ch_metadata_bt2_spikein | view

    /*
     *  SUBWORKFLOW: Filter reads based some standard measures
     *  - Unmapped reads 0x004
     *  - Mate unmapped 0x0008
     *  - Multi-mapped reads
     *  - Filter out reads aligned to blacklist regions
     *  - Filter out reads below a threshold q score
     *  - Filter out mitochondrial reads (if required)
     */
    if (params.run_read_filter) {
        FILTER_READS (
            ch_samtools_bam,
            PREPARE_GENOME.out.allowed_regions.collect{it[1]}.ifEmpty([]),
            PREPARE_GENOME.out.fasta
        )
        ch_samtools_bam      = FILTER_READS.out.bam
        ch_samtools_bai      = FILTER_READS.out.bai
        ch_samtools_stats    = FILTER_READS.out.stats
        ch_samtools_flagstat = FILTER_READS.out.flagstat
        ch_samtools_idxstats = FILTER_READS.out.idxstats
        ch_software_versions = ch_software_versions.mix(FILTER_READS.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * MODULE: Run preseq on BAM files before de-duplication
    */
    ch_preseq_output = Channel.empty()
    if (params.run_preseq) {
        PRESEQ_LCEXTRAP (
            ch_samtools_bam
        )
        ch_preseq_output = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

    /*
     * SUBWORKFLOW: Mark duplicates on all samples
     */
    ch_markduplicates_metrics = Channel.empty()
    if (params.run_mark_dups) {
        MARK_DUPLICATES_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            true,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam           = MARK_DUPLICATES_PICARD.out.bam
        ch_samtools_bai           = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_metrics = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Remove duplicates - default is on IgG controls only
     */
    if (params.run_remove_dups) {
        DEDUPLICATE_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            params.dedup_target_reads,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam      = DEDUPLICATE_PICARD.out.bam
        ch_samtools_bai      = DEDUPLICATE_PICARD.out.bai
        ch_samtools_stats    = DEDUPLICATE_PICARD.out.stats
        ch_samtools_flagstat = DEDUPLICATE_PICARD.out.flagstat
        ch_samtools_idxstats = DEDUPLICATE_PICARD.out.idxstats
        ch_software_versions = ch_software_versions.mix(DEDUPLICATE_PICARD.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bai | view

    /*
    * SUBWORKFLOW: extract duplication stats from picard report
    */
    ch_metadata_picard_duplicates = Channel.empty()
    if (params.run_mark_dups) {
        EXTRACT_PICARD_DUP_META (
            ch_markduplicates_metrics,
            ch_dummy_file.collect(),
            false
        )
        ch_metadata_picard_duplicates = EXTRACT_PICARD_DUP_META.out.metadata
        ch_software_versions          = ch_software_versions.mix(EXTRACT_PICARD_DUP_META.out.versions)
    }
    //ch_metadata_picard_duplicates | view


    /*
     * SUBWORKFLOW: Remove linear amplification duplicates - default is false
     */
    ch_linear_metrics         = Channel.empty()
    ch_linear_duplication_mqc = Channel.empty()
    if (params.run_remove_linear_dups) {
        DEDUPLICATE_LINEAR (
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect(),
            params.dedup_target_reads,
            ch_linear_duplication_header_multiqc
        )
        ch_samtools_bam           = DEDUPLICATE_LINEAR.out.bam
        ch_samtools_bai           = DEDUPLICATE_LINEAR.out.bai
        ch_samtools_stats         = DEDUPLICATE_LINEAR.out.stats
        ch_samtools_flagstat      = DEDUPLICATE_LINEAR.out.flagstat
        ch_samtools_idxstats      = DEDUPLICATE_LINEAR.out.idxstats
        ch_linear_metrics         = DEDUPLICATE_LINEAR.out.metrics
        ch_linear_duplication_mqc = DEDUPLICATE_LINEAR.out.linear_metrics_mqc
        ch_software_versions      = ch_software_versions.mix(DEDUPLICATE_LINEAR.out.versions)
    }

    ch_bedgraph               = Channel.empty()
    ch_bigwig                 = Channel.empty()
    ch_seacr_peaks            = Channel.empty()
    ch_macs2_peaks            = Channel.empty()
    ch_peaks_primary          = Channel.empty()
    ch_peaks_secondary        = Channel.empty()
    ch_peaks_summits          = Channel.empty()
    ch_consensus_peaks        = Channel.empty()
    ch_consensus_peaks_unfilt = Channel.empty()
    if(params.run_peak_calling) {
        /*
        * SUBWORKFLOW: Convert BAM files to bedgraph/bigwig and apply configured normalisation strategy
        */
        PREPARE_PEAKCALLING(
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.chrom_sizes.collect(),
            ch_dummy_file,
            params.normalisation_mode,
            ch_metadata_bt2_spikein
        )
        ch_bedgraph          = PREPARE_PEAKCALLING.out.bedgraph
        ch_bigwig            = PREPARE_PEAKCALLING.out.bigwig
        ch_software_versions = ch_software_versions.mix(PREPARE_PEAKCALLING.out.versions)

        /*
         * CHANNEL: Separate bedgraphs into target/control
         */
        ch_bedgraph.filter { it -> it[0].is_control == false }
        .set { ch_bedgraph_target }
        ch_bedgraph.filter { it -> it[0].is_control == true }
        .set { ch_bedgraph_control }
        //ch_bedgraph_target | view
        //ch_bedgraph_control | view

        /*
        * CHANNEL: Separate bams into target/control
        */
        ch_samtools_bam.filter { it -> it[0].is_control == false }
        .set { ch_bam_target }
        ch_samtools_bam.filter { it -> it[0].is_control == true }
        .set { ch_bam_control }
        //ch_bam_target | view
        //ch_bam_control | view

        if(params.use_control) {
            /*
            * MODULE: Call peaks using SEACR with IgG control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bedgraph_control.map{ row -> [row[0].control_group + "_" + row[0].replicate, row] }
                .cross( ch_bedgraph_target.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bedgraph_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

                SEACR_CALLPEAK_IGG (
                    ch_bedgraph_paired,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks       = SEACR_CALLPEAK_IGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_IGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_CALLPEAK_IGG.out.bed | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bam_control.map{ row -> [row[0].control_group + "_" + row[0].replicate, row] }
                .cross( ch_bam_target.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bam_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BAM, CONTROL_BAM]
                //ch_bam_paired | view

                MACS2_CALLPEAK_IGG (
                    ch_bam_paired,
                    params.macs_gsize
                )
                ch_macs2_peaks       = MACS2_CALLPEAK_IGG.out.peak
                ch_peaks_summits     = MACS2_CALLPEAK_IGG.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_IGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_CALLPEAK_IGG.out.peak | view
            }
        }
        else {
            /*
            * MODULE: Call peaks without IgG Control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bedgraph_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_bedgraph_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
                // ch_bedgraph_target_fctrl | view

                SEACR_CALLPEAK_NOIGG (
                    ch_bedgraph_target_fctrl,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks       = SEACR_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_NO_IGG.out.bed | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bam_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_samtools_bam_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]
                //ch_samtools_bam_target_fctrl | view

                MACS2_CALLPEAK_NOIGG (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks       = MACS2_CALLPEAK_NOIGG.out.peak
                ch_peaks_summits     = MACS2_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                // MACS2_CALLPEAK_NOIGG.out.peak | view
            }
        }

        if ("macs2" in params.callers) {
            /*
            * MODULE: Convert narrow or broad peak to bed
            */
            PEAK_TO_BED ( ch_macs2_peaks )
            ch_macs2_peaks       = PEAK_TO_BED.out.file
            ch_software_versions = ch_software_versions.mix(PEAK_TO_BED.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //PEAK_TO_BED.out.file | view
        }

        // Identify the primary peak data stream for downstream analysis
        if(callers[0] == 'seacr') {
            ch_peaks_primary   = ch_seacr_peaks
            ch_peaks_secondary = ch_macs2_peaks
        }
        if(callers[0] == 'macs2') {
            ch_peaks_primary   = ch_macs2_peaks
            ch_peaks_secondary = ch_seacr_peaks
        }

        if(callers[0] == 'seacr') {
            /*
            * MODULE: Extract summits from seacr peak beds
            */
            AWK_EXTRACT_SUMMITS (
                ch_peaks_primary
            )
            ch_peaks_summits     = AWK_EXTRACT_SUMMITS.out.file
            ch_software_versions = ch_software_versions.mix(AWK_EXTRACT_SUMMITS.out.versions)
            //AWK_EXTRACT_SUMMITS.out.file | view
        }

        /*
        * MODULE: Add sample identifier column to peak beds
        */
        AWK_NAME_PEAK_BED (
            ch_peaks_primary
        )
        ch_software_versions = ch_software_versions.mix(AWK_NAME_PEAK_BED.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //AWK_NAME_PEAK_BED.out.file | view

        if(params.run_consensus_all) {
            /*
            * CHANNEL: Group all samples, filter where the number in the group is > 1
            */
            AWK_NAME_PEAK_BED.out.file
            .map { row -> [ 1, row[1] ] }
            .groupTuple(by: [0])
            .map { row ->
                def new_meta = [:]
                new_meta.put( "id", "all_samples" )
                [ new_meta, row[1].flatten() ]
            }
            .map { row ->
                [ row[0], row[1], row[1].size() ]
            }
            .filter { row -> row[2] > 1 }
            .map { row ->
                [ row[0], row[1] ]
            }
            .set { ch_peaks_bed_all }
            // EXAMPLE CHANNEL STRUCT: [[id: all_samples], [BED1, BED2, BEDn...], count]
            //ch_peaks_bed_all | view

            /*
            * SUBWORKFLOW: Construct group consensus peaks
            */
            CONSENSUS_PEAKS_ALL (
                ch_peaks_bed_all
            )
            ch_consensus_peaks        = CONSENSUS_PEAKS_ALL.out.filtered_bed
            ch_consensus_peaks_unfilt = CONSENSUS_PEAKS_ALL.out.merged_bed
            ch_software_versions      = ch_software_versions.mix(CONSENSUS_PEAKS_ALL.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS_ALL.out.bed | view
        } else {
            /*
            * CHANNEL: Group samples based on group name
            */
            AWK_NAME_PEAK_BED.out.file
            .map { row -> [ row[0].group, row[1] ] }
            .groupTuple(by: [0])
            .map { row -> [ [id: row[0]], row[1].flatten() ] }
            .set { ch_peaks_bed_group }
            // EXAMPLE CHANNEL STRUCT: [[id: <GROUP>], [BED1, BED2, BEDn...], count]
            //ch_peaks_bed_group | view

            /*
            * SUBWORKFLOW: Construct group consensus peaks
            * where there is more than 1 replicate in a group
            */
            CONSENSUS_PEAKS (
                ch_peaks_bed_group
            )
            ch_consensus_peaks        = CONSENSUS_PEAKS.out.filtered_bed
            ch_consensus_peaks_unfilt = CONSENSUS_PEAKS.out.merged_bed
            ch_software_versions      = ch_software_versions.mix(CONSENSUS_PEAKS.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS.out.bed | view
        }
    }

    ch_dt_corrmatrix              = Channel.empty()
    ch_dt_pcadata                 = Channel.empty()
    ch_dt_fpmatrix                = Channel.empty()
    ch_peakqc_frip_mqc            = Channel.empty()
    ch_peakqc_count_mqc           = Channel.empty()
    ch_peakqc_count_consensus_mqc = Channel.empty()
    ch_peakqc_reprod_perc_mqc     = Channel.empty()
    ch_frag_len_hist_mqc          = Channel.empty()
    if(params.run_reporting) {
        if(params.run_igv) {
            /*
            * MODULE: Create igv session
            */
            IGV_SESSION (
                PREPARE_GENOME.out.fasta.map {it[1]},
                PREPARE_GENOME.out.fasta_index.map {it[1]},
                PREPARE_GENOME.out.bed_index,
                //PREPARE_GENOME.out.gtf.collect(),
                ch_peaks_primary.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                ch_peaks_secondary.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                ch_bigwig.collect{it[1]}.ifEmpty([]),
                params.igv_sort_by_groups
            )
            //ch_software_versions = ch_software_versions.mix(IGV_SESSION.out.versions)
        }

        if (params.run_deeptools_heatmaps && params.run_peak_calling) {
            /*
            * CHANNEL: Remove IgG from bigwig channel
            */
            ch_bigwig.filter { it[0].is_control == false }
            .set { ch_bigwig_no_igg }
            // ch_bigwig_no_igg | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Genes
            */
            DEEPTOOLS_COMPUTEMATRIX_GENE (
                ch_bigwig_no_igg,
                PREPARE_GENOME.out.bed.collect()
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX_GENE.out.versions)

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_GENE (
                DEEPTOOLS_COMPUTEMATRIX_GENE.out.matrix
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTHEATMAP_GENE.out.versions)

            /*
            * CHANNEL: Structure output for join on id
            */
            ch_peaks_summits
            .map { row -> [row[0].id, row ].flatten()}
            .set { ch_peaks_summits_id }
            //ch_peaks_bed_id | view

            /*
            * CHANNEL: Join beds and bigwigs on id
            */
            ch_bigwig_no_igg
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_peaks_summits_id )
            .filter ( it -> it[-1].size() > 1)
            .set { ch_dt_bigwig_summits }
            //ch_dt_peaks | view

            ch_dt_bigwig_summits
            .map { row -> row[1,2] }
            .set { ch_ordered_bigwig }
            //ch_ordered_bigwig | view

            ch_dt_bigwig_summits
            .map { row -> row[-1] }
            .set { ch_ordered_peaks_max }
            //ch_ordered_peaks_max | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Peaks
            */

            DEEPTOOLS_COMPUTEMATRIX_PEAKS (
                ch_ordered_bigwig,
                ch_ordered_peaks_max
            )

            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.versions)
            //EXAMPLE CHANNEL STRUCT: [[META], MATRIX]
            //DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix | view

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_PEAKS (
                DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTHEATMAP_PEAKS.out.versions)

            if(params.dt_calc_all_matrix) {
                /*
                * MODULE: Run calc gene matrix for all samples
                */
                DEEPTOOLS_COMPUTEMATRIX_GENE_ALL (
                    ch_bigwig_no_igg.map{it[1]}.toSortedList().map{ [[id:'all_genes'], it]},
                    PREPARE_GENOME.out.bed.toSortedList()
                )

                /*
                * MODULE: Calculate DeepTools heatmap for all samples
                */
                DEEPTOOLS_PLOTHEATMAP_GENE_ALL (
                    DEEPTOOLS_COMPUTEMATRIX_GENE_ALL.out.matrix
                )
            }
        }

        if(params.run_deeptools_qc) {
            /*
            * SUBWORKFLOW: Run suite of deeptools QC on bam files
            */
            DEEPTOOLS_QC (
                ch_samtools_bam,
                ch_samtools_bai,
                params.dt_qc_corr_method
            )
            ch_dt_corrmatrix     = DEEPTOOLS_QC.out.correlation_matrix
            ch_dt_pcadata        = DEEPTOOLS_QC.out.pca_data
            ch_dt_fpmatrix       = DEEPTOOLS_QC.out.fingerprint_matrix
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_QC.out.versions)
        }

        /*
        * CHANNEL: Filter bais for target only
        */
        ch_samtools_bai.filter { it -> it[0].is_control == false }
        .set { ch_bai_target }
        //ch_bai_target | view

        if (params.run_peak_qc && params.run_peak_calling) {
            /*
            * CHANNEL: Filter flagstat for target only
            */
            ch_samtools_flagstat.filter { it -> it[0].is_control == false }
            .set { ch_flagstat_target }
            //ch_flagstat_target | view

            /*
            * SUBWORKFLOW: Extract fragments from bam files for fragment-based FRiP score
            */
            EXTRACT_FRAGMENTS (
                ch_bam_target
            )

            /*
            * SUBWORKFLOW: Run suite of peak QC on peaks
            */
            PEAK_QC(
                ch_peaks_primary,
                AWK_NAME_PEAK_BED.out.file,
                ch_consensus_peaks,
                ch_consensus_peaks_unfilt,
                EXTRACT_FRAGMENTS.out.bed,
                ch_flagstat_target,
                params.min_frip_overlap,
                ch_frip_score_header_multiqc,
                ch_peak_counts_header_multiqc,
                ch_peak_counts_consensus_header_multiqc,
                ch_peak_reprod_header_multiqc
            )
            ch_peakqc_frip_mqc             = PEAK_QC.out.primary_frip_mqc
            ch_peakqc_count_mqc            = PEAK_QC.out.primary_count_mqc
            ch_peakqc_count_consensus_mqc  = PEAK_QC.out.consensus_count_mqc
            ch_peakqc_reprod_perc_mqc      = PEAK_QC.out.reprod_perc_mqc
            ch_software_versions           = ch_software_versions.mix(PEAK_QC.out.versions)
        }
        //ch_peakqc_reprod_perc_mqc | view

        /*
        * CHANNEL: Combine bam and bai files on id
        */

        ch_bam_target.map { row -> [row[0].id, row ].flatten()}
        .join ( ch_bai_target.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bam_bai | view

        /*
        * MODULE: Calculate fragment lengths
        */
        SAMTOOLS_CUSTOMVIEW (
            ch_bam_bai
        )
        ch_software_versions = ch_software_versions.mix(SAMTOOLS_CUSTOMVIEW.out.versions)
        //SAMTOOLS_CUSTOMVIEW.out.tsv | view

        /*
        * CHANNEL: Prepare data for generate reports
        */
        // Make sure files are always in order for resume
        ch_frag_len = SAMTOOLS_CUSTOMVIEW.out.tsv
        .toSortedList { row -> row[0].id }
        .map {
            list ->
            def output = []
            list.each{ v -> output.add(v[1]) }
            output
        }
        //ch_frag_len | view

        /*
        * MODULE: Calculate fragment length histogram for mqc
        */
        FRAG_LEN_HIST(
            ch_frag_len,
            ch_frag_len_header_multiqc
        )
        ch_frag_len_hist_mqc = FRAG_LEN_HIST.out.frag_len_mqc
        ch_software_versions = ch_software_versions.mix(FRAG_LEN_HIST.out.versions)
    }
    //ch_frag_len_hist_mqc | view


    if (params.run_multiqc) {
        workflow_summary    = WorkflowCutandrun.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        /*
        * MODULE: Collect software versions used in pipeline
        */
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_software_versions.unique().collectFile()
        )

        /*
        * MODULE: Multiqc
        */
        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect(),
            ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yml"),
            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([]),
            ch_preseq_output.collect{it[1]}.ifEmpty([]),
            ch_dt_corrmatrix.collect{it[1]}.ifEmpty([]),
            ch_dt_pcadata.collect{it[1]}.ifEmpty([]),
            ch_dt_fpmatrix.collect{it[1]}.ifEmpty([]),
            ch_peakqc_count_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_frip_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_count_consensus_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_reprod_perc_mqc.collect().ifEmpty([]),
            ch_frag_len_hist_mqc.collect().ifEmpty([]),
            ch_linear_duplication_mqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
