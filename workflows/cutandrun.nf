/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Create summary input parameters map for reporting
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters in specialised library
WorkflowCutandrun.initialise(params, log)

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.blacklist,
    params.bowtie2,
    params.fasta,
    params.gtf,
    params.input,
    params.spikein_bowtie2,
    params.spikein_fasta
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

ch_blacklist = Channel.empty()
if (params.blacklist) {
    ch_blacklist = file(params.blacklist)
}
else {
    ch_blacklist = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
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
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_frag_len_header_multiqc = file("$projectDir/assets/multiqc/frag_len_header.txt", checkIfExists: true)

/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

// Init aligners
def prepare_tool_indices = ["bowtie2"]

// Check normalisation mode params
def norm_mode_list = ["Spikein", "RPKM", "CPM", "BPM", "None" ]
if (!(params.normalisation_mode in norm_mode_list)) {
    exit 1, "Invalid normalisation mode option: ${params.normalisation_mode}. Valid options: ${norm_mode_list.join(', ')}"
}

// Check peak caller params
def caller_list = ['seacr', 'macs2']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : ['seacr']
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

// Check consensus peak mode params
def conseneus_mode_list = ['group', 'all']
if (!(params.consensus_peak_mode in conseneus_mode_list)) {
    exit 1, "Invalid conseneus mode option: ${params.consensus_peak_mode}. Valid options: ${conseneus_mode_list.join(', ')}"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { INPUT_CHECK                     } from "../subworkflows/local/input_check"
include { AWK as AWK_NAME_PEAK_BED        } from "../modules/local/linux/awk"
include { AWK as AWK_FRAG_BIN             } from "../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW             } from "../modules/local/samtools_custom_view"
include { IGV_SESSION                     } from "../modules/local/python/igv_session"
include { AWK as AWK_EDIT_PEAK_BED        } from "../modules/local/linux/awk"
include { CALCULATE_FRIP                  } from "../modules/local/modules/calculate_frip/main"
include { CUT as CUT_CALC_REPROD          } from "../modules/local/linux/cut"
include { CALCULATE_PEAK_REPROD           } from "../modules/local/modules/calculate_peak_reprod/main"
include { EXPORT_META                     } from "../modules/local/export_meta"
include { EXPORT_META as EXPORT_META_CTRL } from "../modules/local/export_meta"
include { GENERATE_REPORTS                } from "../modules/local/modules/generate_reports/main"
include { MULTIQC                         } from "../modules/local/multiqc"

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
include { PREPARE_GENOME                                 } from "../subworkflows/local/prepare_genome"
include { FASTQC_TRIMGALORE                              } from "../subworkflows/local/fastqc_trimgalore"
include { ALIGN_BOWTIE2                                  } from "../subworkflows/local/align_bowtie2"
include { ANNOTATE_META_AWK as ANNOTATE_BT2_META         } from "../subworkflows/local/annotate_meta_awk"
include { ANNOTATE_META_AWK as ANNOTATE_BT2_SPIKEIN_META } from "../subworkflows/local/annotate_meta_awk"
include { ANNOTATE_META_AWK as ANNOTATE_DEDUP_META       } from "../subworkflows/local/annotate_meta_awk"
include { CONSENSUS_PEAKS                                } from "../subworkflows/local/consensus_peaks"
include { CONSENSUS_PEAKS as CONSENSUS_PEAKS_ALL         } from "../subworkflows/local/consensus_peaks"
include { CALCULATE_FRAGMENTS                            } from "../subworkflows/local/calculate_fragments"
include { ANNOTATE_META_CSV as ANNOTATE_FRIP_META        } from "../subworkflows/local/annotate_meta_csv"
include { ANNOTATE_META_CSV as ANNOTATE_PEAK_REPRO_META  } from "../subworkflows/local/annotate_meta_csv"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CAT_FASTQ                                                } from "../modules/nf-core/modules/cat/fastq/main"
include { SEACR_CALLPEAK                                           } from "../modules/nf-core/modules/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                   } from "../modules/nf-core/modules/seacr/callpeak/main"
include { MACS2_CALLPEAK                                           } from "../modules/nf-core/modules/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG                   } from "../modules/nf-core/modules/macs2/callpeak/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE  } from "../modules/nf-core/modules/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS } from "../modules/nf-core/modules/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE      } from "../modules/nf-core/modules/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS     } from "../modules/nf-core/modules/deeptools/plotheatmap/main"
include { BEDTOOLS_INTERSECT                                       } from "../modules/nf-core/modules/bedtools/intersect/main.nf"
include { CUSTOM_DUMPSOFTWAREVERSIONS                              } from "../modules/local/modules/custom/dumpsoftwareversions/main"

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { MARK_DUPLICATES_PICARD                       } from "../subworkflows/nf-core/mark_duplicates_picard"
include { MARK_DUPLICATES_PICARD as DEDUPLICATE_PICARD } from "../subworkflows/nf-core/mark_duplicates_picard"
include { SAMTOOLS_VIEW_SORT_STATS                     } from "../subworkflows/nf-core/samtools_view_sort_stats"
include { PREPARE_PEAKCALLING                          } from "../subworkflows/nf-core/prepare_peakcalling"

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CUTANDRUN {

    // Init
    ch_software_versions = Channel.empty()
    ch_frag_len_multiqc  = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    if(params.run_genome_prep) {
        PREPARE_GENOME (
            prepare_tool_indices
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
                PREPARE_GENOME.out.bowtie2_spikein_index
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
    // ch_samtools_bam | view

    /*
     *  SUBWORKFLOW: Filter reads based on quality metrics
     *  http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
     */
    if (params.run_q_filter) {
        SAMTOOLS_VIEW_SORT_STATS (
            ch_samtools_bam
        )
        ch_samtools_bam      = SAMTOOLS_VIEW_SORT_STATS.out.bam
        ch_samtools_bai      = SAMTOOLS_VIEW_SORT_STATS.out.bai
        ch_samtools_stats    = SAMTOOLS_VIEW_SORT_STATS.out.stats
        ch_samtools_flagstat = SAMTOOLS_VIEW_SORT_STATS.out.flagstat
        ch_samtools_idxstats = SAMTOOLS_VIEW_SORT_STATS.out.idxstats
        ch_software_versions = ch_software_versions.mix(SAMTOOLS_VIEW_SORT_STATS.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Mark duplicates on all samples
     */
    ch_markduplicates_metrics = Channel.empty()
    if (params.run_mark_dups) {
        MARK_DUPLICATES_PICARD (
            ch_samtools_bam,
            true
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
    ch_dedup_multiqc = Channel.empty()
    if (params.run_remove_dups) {
        DEDUPLICATE_PICARD (
            ch_samtools_bam,
            params.dedup_target_reads
        )
        ch_samtools_bam      = DEDUPLICATE_PICARD.out.bam
        ch_samtools_bai      = DEDUPLICATE_PICARD.out.bai
        ch_samtools_stats    = DEDUPLICATE_PICARD.out.stats
        ch_samtools_flagstat = DEDUPLICATE_PICARD.out.flagstat
        ch_samtools_idxstats = DEDUPLICATE_PICARD.out.idxstats
        ch_dedup_multiqc     = DEDUPLICATE_PICARD.out.metrics
        ch_software_versions = ch_software_versions.mix(DEDUPLICATE_PICARD.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Annotate meta-data with aligner stats for target and spike-in
     * the meta-data is annotated additivley so we only need to track the final channel output
     */
    if (params.aligner == "bowtie2" && params.run_alignment) {
        ANNOTATE_BT2_META (
            ch_samtools_bam,
            ch_bowtie2_log,
            ch_bt2_to_csv_awk,
            "",
            "_target",
            true
        )
        ch_software_versions = ch_software_versions.mix(ANNOTATE_BT2_META.out.versions)

        ANNOTATE_BT2_SPIKEIN_META (
            ANNOTATE_BT2_META.out.output,
            ch_bowtie2_spikein_log,
            ch_bt2_to_csv_awk,
            "",
            "_spikein",
            true
        )
        ch_samtools_bam = ANNOTATE_BT2_SPIKEIN_META.out.output
    }
    // META-DATA example state:
    //[[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false,
    // bt2_total_reads_spikein:9616, bt2_align1_spikein:1, bt2_align_gt1_spikein:0, bt2_non_aligned_spikein:9615, bt2_total_aligned_spikein:1,
    // bt2_total_reads_target:9616, bt2_align1_target:315, bt2_align_gt1_target:449, bt2_non_aligned_target:8852, bt2_total_aligned_target:764], BAM]
    //ch_samtools_bam | view
    //EXPORT_META ( ch_annotated_meta.collect{ it[0] } )

    /*
    * SUBWORKFLOW: Annotate meta-data with duplication stats
    */
    if (params.run_mark_dups) {
        ANNOTATE_DEDUP_META(
            ch_samtools_bam,
            ch_markduplicates_metrics,
            ch_dummy_file.collect(),
            "dedup_",
            "",
            false
        )
        ch_samtools_bam      = ANNOTATE_DEDUP_META.out.output
        ch_software_versions = ch_software_versions.mix(ANNOTATE_DEDUP_META.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[META + dedup_library:unknown library, dedup_unpaired_reads_examined:0, dedup_read_pairs_examined:350, dedup_secondary_or_supplementary_rds:0,
    // dedup_unmapped_reads:0, dedup_unpaired_read_duplicates:0, dedup_read_pair_duplicates:0, dedup_read_pair_optical_duplicates:0, dedup_percent_duplication:0,
    // dedup_estimated_library_size:], BAM]
    //ch_samtools_bam | view

    ch_bedgraph = Channel.empty()
    ch_bigwig   = Channel.empty()
    if(params.run_peak_calling) {
        /*
        * SUBWORKFLOW: Convert bam files to bedgraph/bigwig and apply configured normalisation strategy
        */
        PREPARE_PEAKCALLING(
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.chrom_sizes,
            ch_dummy_file,
            params.normalisation_mode,
        )
        ch_samtools_bam      = PREPARE_PEAKCALLING.out.bam
        ch_bedgraph          = PREPARE_PEAKCALLING.out.bedgraph
        ch_bigwig            = PREPARE_PEAKCALLING.out.bigwig
        ch_software_versions = ch_software_versions.mix(ANNOTATE_DEDUP_META.out.versions)

        /*
         * CHANNEL: Separate bedgraphs into target/control
         */
        ch_bedgraph.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_bedgraph_split }
        //ch_bedgraph_split.target | view
        //ch_bedgraph_split.control | view

        ch_seacr_bed = Channel.empty()
        ch_macs2_bed = Channel.empty()
        ch_peaks_bed = Channel.empty()

        if(params.use_control) {
            /*
             * MODULE: Call peaks using SEACR with IgG control
             */
            if('seacr' in callers) {
                /*
                * CHANNEL: Pull control groups
                */
                ch_bedgraph_split.target.map{
                    row -> [row[0].control_group, row]
                }
                .set { ch_bg_target_ctrlgrp }
                //ch_bg_target_ctrlgrp | view

                ch_bedgraph_split.control.map{
                    row -> [row[0].control_group, row]
                }
                .set { ch_bg_control_ctrlgrp }
                //ch_bg_control_ctrlgrp | view

                /*
                * CHANNEL: Create target/control pairings
                */
                // Create pairs of controls (IgG) with target samples if they are supplied
                ch_bg_control_ctrlgrp.cross(ch_bg_target_ctrlgrp).map {
                    row -> [row[1][1][0], row[1][1][1], row[0][1][1]]
                }
                .set{ ch_bedgraph_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]
                //ch_bedgraph_paired | view

                SEACR_CALLPEAK (
                    ch_bedgraph_paired,
                    params.peak_threshold
                )
                ch_seacr_bed         = SEACR_CALLPEAK.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_CALLPEAK.out.bed | view
            }

            if('macs2' in callers) {
                ch_samtools_bam.branch{ it ->
                    target:  it[0].is_control == false
                    control: it[0].is_control == true
                }
                .set { ch_samtools_bam_split }
                // ch_samtools_bam_split.target | view

                /*
                * CHANNEL: Pull control groups
                */
                ch_samtools_bam_split.target.map{
                    row -> [row[0].control_group, row]
                }
                .set { ch_bam_target_ctrlgrp }
                //ch_bam_target_ctrlgrp | view

                ch_samtools_bam_split.control.map{
                    row -> [row[0].control_group, row]
                }
                .set { ch_bam_control_ctrlgrp }
                // ch_bam_control_ctrlgrp | view

                /*
                * CHANNEL: Create target/control pairings
                */
                // Create pairs of controls (IgG) with target samples if they are supplied
                ch_bam_control_ctrlgrp.cross(ch_bam_target_ctrlgrp).map{
                    row -> [row[1][1][0], row[1][1][1], row[0][1][1]]
                }
                .set{ch_bam_paired}
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BAM, CONTROL_BAM]
                // ch_bam_paired | view

                MACS2_CALLPEAK (
                    ch_bam_paired,
                    params.macs2_gsize
                )
                ch_macs2_bed         = MACS2_CALLPEAK.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_CALLPEAK.out.bed | view
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
                ch_bedgraph_split.target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_bedgraph_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
                // ch_bedgraph_target_fctrl | view

                SEACR_CALLPEAK_NOIGG (
                    ch_bedgraph_target_fctrl,
                    params.peak_threshold
                )
                ch_seacr_bed         = SEACR_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_NO_IGG.out.bed | view
            }

            if('macs2' in callers) {
                ch_samtools_bam.branch{ it ->
                    target:  it[0].is_control == false
                    control: it[0].is_control == true
                }
                .set { ch_samtools_bam_split }
                // ch_samtools_bam_split.target | view

                /*
                * CHANNEL: Add fake control channel
                */
                ch_samtools_bam_split.target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_samtools_bam_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]
                //ch_samtools_bam_target_fctrl | view

                MACS2_CALLPEAK_NOIGG (
                    ch_samtools_bam_target_fctrl,
                    params.macs2_gsize
                )
                ch_macs2_bed         = MACS2_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                // MACS2_CALLPEAK_NOIGG.out.bed | view
            }
        }

        // Store output of primary peakcaller in the output channel
        if(callers[0] == 'seacr') {
                ch_peaks_bed = ch_seacr_bed
        }
        if(callers[0] == 'macs2') {
            ch_peaks_bed = ch_macs2_bed
        }

        /*
        * MODULE: Add sample identifier column to peak beds
        */
        AWK_NAME_PEAK_BED (
            ch_peaks_bed
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
                new_meta = [:]
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
                ch_peaks_bed_all,
                params.skip_upset_plots
            )
            ch_software_versions = ch_software_versions.mix(CONSENSUS_PEAKS_ALL.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS_ALL.out.bed | view

        } else {
            /*
            * CHANNEL: Group samples based on group
            */
            AWK_NAME_PEAK_BED.out.file
            .map { row -> [ row[0].group, row[1] ] }
            .groupTuple(by: [0])
            .map { row ->
                new_meta = [:]
                new_meta.put( "id", row[0] )
                [ new_meta, row[1].flatten() ]
            }
            .map { row ->
                [ row[0], row[1], row[1].size() ]
            }
            .filter { row -> row[2] > 1 }
            .map { row ->
                [ row[0], row[1] ]
            }
            .set { ch_peaks_bed_group }
            // EXAMPLE CHANNEL STRUCT: [[id: <GROUP>], [BED1, BED2, BEDn...], count]
            //ch_peaks_bed_group | view

            /*
            * SUBWORKFLOW: Construct group consensus peaks
            * where there is more than 1 replicate in a group
            */
            CONSENSUS_PEAKS (
                ch_peaks_bed_group,
                params.skip_upset_plots
            )
            ch_software_versions = ch_software_versions.mix(CONSENSUS_PEAKS.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS.out.bed | view
        }

        /*
        * SUBWORKFLOW: Calculate fragment bed from bams
        * - Filter for mapped reads
        * - Convert to bed file
        * - Keep the read pairs that are on the same chromosome and fragment length less than 1000bp
        * - Only extract the fragment related columns using cut
        */
        CALCULATE_FRAGMENTS (
            ch_samtools_bam
        )
        ch_software_versions = ch_software_versions.mix(CALCULATE_FRAGMENTS.out.versions)
        //EXAMPLE CHANNEL STRUCT: NO CHANGE
        //CALCULATE_FRAGMENTS.out.bed | view

        /*
        * MODULE: Bin the fragments into 500bp bins ready for downstream reporting
        */
        AWK_FRAG_BIN(
            CALCULATE_FRAGMENTS.out.bed
        )
        //AWK_FRAG_BIN.out.file | view


        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_samtools_bam.map { row -> [row[0].id, row ].flatten()}
        .join ( ch_samtools_bai.map { row -> [row[0].id, row ].flatten()} )
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
        // ch_software_versions = ch_software_versions.mix(SAMTOOLS_CUSTOMVIEW.out.versions)
        //SAMTOOLS_CUSTOMVIEW.out.tsv | view
    }

    if(params.run_reporting) {
        if(params.run_igv) {
            /*
            * MODULE: Create igv session
            */
            IGV_SESSION (
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.gtf,
                ch_peaks_bed.collect{it[1]}.ifEmpty([]),
                ch_bigwig.collect{it[1]}.ifEmpty([])
            )
            //ch_software_versions = ch_software_versions.mix(IGV_SESSION.out.versions)
        }

        if (params.run_deep_tools && params.run_peak_calling) {
            /*
            * MODULE: Extract max signal from peak beds
            */
            AWK_EDIT_PEAK_BED (
                ch_peaks_bed
            )
            ch_software_versions = ch_software_versions.mix(AWK_EDIT_PEAK_BED.out.versions)
            //AWK_EDIT_PEAK_BED.out.file | view

            /*
            * CHANNEL: Structure output for join on id
            */
            AWK_EDIT_PEAK_BED.out.file
                .map { row -> [row[0].id, row ].flatten()}
                .set { ch_peaks_bed_id }
            //ch_peaks_bed_id | view

            /*
            * CHANNEL: Remove IgG from bigwig channel
            */
            ch_bigwig
                .filter { it[0].is_control == false }
                .set { ch_bigwig_no_igg }
            //ch_bigwig_no_igg | view

            /*
            * CHANNEL: Join beds and bigwigs on id
            */
            ch_bigwig_no_igg
                .map { row -> [row[0].id, row ].flatten()}
                .join ( ch_peaks_bed_id )
                .set { ch_dt_peaks }
            //ch_dt_peaks | view

            ch_dt_peaks
                .map { row -> row[1,2] }
                .set { ch_ordered_bigwig }
            //ch_ordered_bigwig | view

            ch_dt_peaks
                .map { row -> row[-1] }
                .set { ch_ordered_peaks_max }
            //ch_ordered_peaks_max | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Genes
            */
            DEEPTOOLS_COMPUTEMATRIX_GENE (
                ch_bigwig_no_igg,
                PREPARE_GENOME.out.bed
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX_GENE.out.versions)

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_GENE (
                DEEPTOOLS_COMPUTEMATRIX_GENE.out.matrix
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTHEATMAP_GENE.out.versions)

            // Run if not empty file size > 1 byte
            ch_ordered_peaks_max
                .filter { it -> it.size() > 1}
                .set { ch_ordered_peaks_max_notempty }
            //ch_ordered_peaks_max_notempty | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Peaks
            */
            DEEPTOOLS_COMPUTEMATRIX_PEAKS (
                ch_ordered_bigwig,
                ch_ordered_peaks_max_notempty
            )
            //EXAMPLE CHANNEL STRUCT: [[META], MATRIX]
            //DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix | view

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_PEAKS (
                DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
            )
        }

        /*
        * CHANNEL: Join bams and beds on id
        */
        ch_samtools_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_samtools_bai.map { row -> [row[0].id, row ].flatten()} )
            .join ( ch_peaks_bed.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4], row[6]] }
            .set { ch_bam_bai_bed }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, BED]
        //ch_bam_bai_bed | view

        ch_samtools_bam_ctrl = ch_samtools_bam
        if(!params.skip_frip) {
            /*
            * MODULE: Calculate Frip scores for samples
            */
            CALCULATE_FRIP (
                ch_bam_bai_bed
            )
            ch_software_versions = ch_software_versions.mix(CALCULATE_FRIP.out.versions)

            /*
            * SUBWORKFLOW: Annotate meta-data with frip stats
            */
            ANNOTATE_FRIP_META (
                ch_samtools_bam,
                CALCULATE_FRIP.out.frips,
                "",
                ""
            )
            ch_samtools_bam = ANNOTATE_FRIP_META.out.output
            //ch_samtools_bam | view
        }
       
        /*
        * MODULE: Trim unwanted columns for downstream reporting
        */
        CUT_CALC_REPROD (
            AWK_NAME_PEAK_BED.out.file
        )

        /*
        * CHANNEL: Group samples based on group
        */
        CUT_CALC_REPROD.out.file
            .map { row -> [ row[0].group, row[1] ] }
            .groupTuple(by: [0])
            .map { row ->
                new_meta = [:]
                new_meta.put( "id", row[0] )
                [ new_meta, row[1].flatten() ]
            }
            .map { row ->
                [ row[0], row[1], row[1].size() ]
            }
            .filter { row -> row[2] > 1 }
            .map { row ->
                [ row[0], row[1] ]
            }
        .set { ch_seacr_bed_group_2 }

        /*
        * CHANNEL: Per group, create a channel per one against all combination
        */
        ch_seacr_bed_group_2
            .flatMap{
                row ->
                new_output = []
                row[1].each{ file ->
                    files_copy = row[1].collect()
                    files_copy.remove(files_copy.indexOf(file))
                    new_output.add([[id: file.name.split("\\.")[0]], file, files_copy])
                }
                new_output
            }
        .set { ch_beds_intersect }
        //EXAMPLE CHANNEL STRUCT: [[META], BED (-a), [BED...n] (-b)]
        //ch_beds_intersect | view

        /*
        * MODULE: Find intra-group overlap
        */
        BEDTOOLS_INTERSECT (
            ch_beds_intersect,
            "bed"
        )
        ch_software_versions = ch_software_versions.mix(BEDTOOLS_INTERSECT.out.versions)
        //EXAMPLE CHANNEL STRUCT: [[META], BED]
        //BEDTOOLS_INTERSECT.out.intersect | view

        /*
        * MODULE: Use overlap to calculate a peak repro %
        */
        CALCULATE_PEAK_REPROD (
            BEDTOOLS_INTERSECT.out.intersect
        )
        ch_software_versions = ch_software_versions.mix(CALCULATE_PEAK_REPROD.out.versions)
        //EXAMPLE CHANNEL STRUCT: [[META], CSV]
        //CALCULATE_PEAK_REPROD.out.csv

        /*
        * SUBWORKFLOW: Annotate meta-data with peak stats
        */
        ANNOTATE_PEAK_REPRO_META (
            ch_samtools_bam,
            CALCULATE_PEAK_REPROD.out.csv,
            "",
            ""
        )
        ch_samtools_bam = ANNOTATE_PEAK_REPRO_META.out.output
        //ch_samtools_bam | view
        //ANNOTATE_PEAK_REPRO_META.out.output | view

        /*
        * MODULE: Export meta-data to csv file
        */
        EXPORT_META (
            ch_samtools_bam.collect{it[0]},
            "meta_table"
        )

        /*
        * MODULE: Export meta-data to csv file
        */
        EXPORT_META_CTRL (
            ch_samtools_bam_ctrl.collect{it[0]},
            "meta_table_ctrl"
        )

        /*
        * MODULE: Generate python reporting using mixture of meta-data and direct file processing
        */
        GENERATE_REPORTS(
            EXPORT_META.out.csv.collect().ifEmpty([]),  // meta-data report stats
            EXPORT_META_CTRL.out.csv,                   // meta-data report stats
            SAMTOOLS_CUSTOMVIEW.out.tsv.collect{it[1]}, // raw fragments
            AWK_FRAG_BIN.out.file.collect{it[1]},       // binned fragments
            ch_peaks_bed.collect{it[1]},                // peak beds
            ch_frag_len_header_multiqc                  // multiqc config header for fragment length distribution plot
        )
        ch_frag_len_multiqc  = GENERATE_REPORTS.out.frag_len_multiqc
        ch_software_versions = ch_software_versions.mix(GENERATE_REPORTS.out.versions)
    }

    /*
    * MODULE: Collect software versions used in pipeline
    */
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile()
    )

    /*
     * MODULE: Multiqc
     */
    if (params.run_multiqc) {
        workflow_summary    = WorkflowCutandrun.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect(),
            ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml"),
            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([]),
            ch_frag_len_multiqc.collect().ifEmpty([])
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
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
