/*
 * Calculate Peak-based metrics and QC
*/

include { PEAK_FRIP                            } from "../../modules/local/peak_frip"
include { PEAK_COUNTS as PRIMARY_PEAK_COUNTS   } from "../../modules/local/peak_counts"
include { PEAK_COUNTS as CONSENSUS_PEAK_COUNTS } from "../../modules/local/peak_counts"

// include { CUT as CUT_CALC_REPROD          } from "../../modules/local/linux/cut"
// include { CALCULATE_PEAK_REPROD           } from "../../modules/local/peak_reprod"
// include { PLOT_CONSENSUS_PEAKS } from '../../modules/local/python/plot_consensus_peaks'

workflow PEAK_QC {
    take:
    peaks                               // channel: [ val(meta), [ bed ] ]
    consensus_peaks                     // channel: [ val(meta), [ bed ] ]
    bam                                 // channel: [ val(meta), [ bam ] ]
    flagstat                            // channel: [ val(meta), [ flagstat ] ]
    min_frip_overlap                    // val
    frip_score_header_multiqc           // file
    peak_count_header_multiqc           // file
    peak_count_consensus_header_multiqc // file

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


    // // Plot consensus peak sets
    // if(!skip_plot) {
    //     ch_merged_bed_sorted = BEDTOOLS_MERGE.out.bed
    //         .toSortedList { row -> row[0].id }
    //         .map { list ->
    //             def output = []
    //             list.each{ v -> output.add(v[1]) }
    //             output
    //         }

    //     PLOT_CONSENSUS_PEAKS ( ch_merged_bed_sorted.ifEmpty([]) )
    //     ch_versions = ch_versions.mix(PLOT_CONSENSUS_PEAKS.out.versions)
    // }

    emit:
    primary_frip_mqc    = PEAK_FRIP.out.frip_mqc              // channel: [ val(meta), [ mqc ] ]
    primary_count_mqc   = PRIMARY_PEAK_COUNTS.out.count_mqc   // channel: [ val(meta), [ mqc ] ]
    consensus_count_mqc = CONSENSUS_PEAK_COUNTS.out.count_mqc // channel: [ val(meta), [ mqc ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
