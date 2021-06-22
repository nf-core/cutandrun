/*
 * Create group consensus peaks
 */

include { SORT } from "../../modules/local/sort" addParams( options: params.sort_options)
include { BEDTOOLS_MERGE } from "../../modules/nf-core/software/bedtools/merge/main" addParams( options: params.bedtools_merge_options )
include { AWK } from "../../modules/local/awk" addParams( options: params.awk_threshold_options)
include { PLOT_CONSENSUS_PEAKS } from "../../modules/local/plot_consensus_peaks" addParams( options: params.plot_peak_options )

workflow CONSENSUS_PEAKS {

    take:
    bed //  channel: [ val(meta), [ bed ] ]

    main:

    // Sort and merge bed files
    SORT ( bed )

    // Merge peaks
    BEDTOOLS_MERGE ( SORT.out.file )

    // Optionally filter peaks on minimum replicate consensus and produce plot
    AWK ( BEDTOOLS_MERGE.out.bed )

    // Plot consensus peak sets
    PLOT_CONSENSUS_PEAKS ( BEDTOOLS_MERGE.out.bed.collect{it[1]}.ifEmpty([]) ) //.collect().ifEmpty([]) )

    emit:
    bed              = BEDTOOLS_MERGE.out.bed       // channel: [ val(meta), [ bed ] ]
    filtered_bed     = AWK.out.file                 // channel: [ val(meta), [ bed ] ]
    bedtools_version = BEDTOOLS_MERGE.out.version   // path: *.version.txt
}
