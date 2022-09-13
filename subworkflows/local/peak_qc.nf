/*
 * Calculate Peak-based metrics and QC
*/

include { PEAK_FRIP                            } from "../../modules/local/peak_frip"
include { PEAK_COUNTS as PRIMARY_PEAK_COUNTS   } from "../../modules/local/peak_counts"
include { PEAK_COUNTS as CONSENSUS_PEAK_COUNTS } from "../../modules/local/peak_counts"
include { CUT as CUT_CALC_REPROD               } from "../../modules/local/linux/cut"
include { BEDTOOLS_INTERSECT                   } from "../../modules/nf-core/modules/bedtools/intersect/main.nf"
include { CALCULATE_PEAK_REPROD                } from "../../modules/local/python/peak_reprod"
include { PLOT_CONSENSUS_PEAKS                 } from '../../modules/local/python/plot_consensus_peaks'

workflow PEAK_QC {
    take:
    peaks                               // channel: [ val(meta), [ bed ] ]
    peaks_with_ids                      // channel: [ val(meta), [ bed ] ]
    consensus_peaks                     // channel: [ val(meta), [ bed ] ]
    consensus_peaks_unfiltered          // channel: [ val(meta), [ bed ] ]
    bam                                 // channel: [ val(meta), [ bam ] ]
    flagstat                            // channel: [ val(meta), [ flagstat ] ]
    min_frip_overlap                    // val
    frip_score_header_multiqc           // file
    peak_count_header_multiqc           // file
    peak_count_consensus_header_multiqc // file
    peak_reprod_header_multiqc          // file

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Calculate frip scores for primary peaks
    */
    PEAK_FRIP(
        peaks,
        bam,
        flagstat,
        frip_score_header_multiqc,
        min_frip_overlap
    )
    ch_versions = ch_versions.mix(PEAK_FRIP.out.versions)
    // PEAK_FRIP.out.frip_mqc | view 

    /*
    * MODULE: Calculate peak counts for primary peaks
    */
    PRIMARY_PEAK_COUNTS(
        peaks,
        peak_count_header_multiqc
    )
    ch_versions = ch_versions.mix(PRIMARY_PEAK_COUNTS.out.versions)
    // PRIMARY_PEAK_COUNTS.out.count_mqc | view

    /*
    * MODULE: Calculate peak counts for consensus peaks
    */
    CONSENSUS_PEAK_COUNTS(
        consensus_peaks,
        peak_count_consensus_header_multiqc
    )
    ch_versions = ch_versions.mix(CONSENSUS_PEAK_COUNTS.out.versions)
    // CONSENSUS_PEAK_COUNTS.out.count_mqc | view

    /*
    * MODULE: Trim unwanted columns for downstream reporting
    */
    CUT_CALC_REPROD (
        peaks_with_ids
    )
    ch_versions = ch_versions.mix(CUT_CALC_REPROD.out.versions)

    /*
    * CHANNEL: Group samples based on group
    */
    CUT_CALC_REPROD.out.file
    .map { row -> [ row[0].group, row[1] ] }
    .groupTuple(by: [0])
    .map { row -> [ [id: row[0]], row[1].flatten() ] }
    .set { ch_peak_bed_group }
    //ch_peak_bed_group | view

    /*
    * CHANNEL: Per group, create a channel per one against all combination
    */
    ch_peak_bed_group.flatMap{
        row ->
        def new_output = []
        row[1].each{ file ->
            def files_copy = row[1].collect()
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
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], BED]
    //BEDTOOLS_INTERSECT.out.intersect | view

    /*
    * MODULE: Use overlap to calculate a peak repro %
    */
    CALCULATE_PEAK_REPROD (
        BEDTOOLS_INTERSECT.out.intersect,
        peak_reprod_header_multiqc
    )
    ch_versions = ch_versions.mix(CALCULATE_PEAK_REPROD.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], TSV]
    //CALCULATE_PEAK_REPROD.out.tsv

    /*
    * CHANNEL: Prep for upset input
    */
    consensus_peaks_unfiltered
    .toSortedList { row -> row[0].id }
    .map { list ->
        def output = []
        list.each{ v -> output.add(v[1]) }
        output
    }
    .set { ch_merged_bed_sorted }

    /*
    * MODULE: Plot upset plots for sample peaks
    */
    PLOT_CONSENSUS_PEAKS ( 
        ch_merged_bed_sorted.ifEmpty([]) 
    )
    ch_versions = ch_versions.mix(PLOT_CONSENSUS_PEAKS.out.versions)

    emit:
    primary_frip_mqc    = PEAK_FRIP.out.frip_mqc              // channel: [ val(meta), [ mqc ] ]
    primary_count_mqc   = PRIMARY_PEAK_COUNTS.out.count_mqc   // channel: [ val(meta), [ mqc ] ]
    consensus_count_mqc = CONSENSUS_PEAK_COUNTS.out.count_mqc // channel: [ val(meta), [ mqc ] ]
    reprod_perc_mqc     = CALCULATE_PEAK_REPROD.out.mqc       // channel: [ val(meta), [ mqc ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
