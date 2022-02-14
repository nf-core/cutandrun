/*
 * Create group consensus peaks
 */


include { SORT                 } from '../../modules/local/linux/sort'
include { BEDTOOLS_MERGE       } from '../../modules/nf-core/modules/bedtools/merge/main'
include { AWK                  } from '../../modules/local/linux/awk'
include { PLOT_CONSENSUS_PEAKS } from '../../modules/local/python/plot_consensus_peaks'

workflow CONSENSUS_PEAKS {

    take:
    bed       //  channel: [ val(meta), [ bed ], count]
    skip_plot // boolean: true/false


    main:
    ch_versions = Channel.empty()

    // Sort bed files
    SORT ( bed )

    // Merge peaks
    BEDTOOLS_MERGE ( SORT.out.file )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    // Filter peaks on minimum replicate consensus
    AWK ( BEDTOOLS_MERGE.out.bed )
    ch_versions = ch_versions.mix(AWK.out.versions)

    // Plot consensus peak sets
    if(!skip_plot) {
        PLOT_CONSENSUS_PEAKS ( BEDTOOLS_MERGE.out.bed.collect{it[1]}.ifEmpty([]) )
        ch_versions = ch_versions.mix(PLOT_CONSENSUS_PEAKS.out.versions)
    }

    emit:
    bed               = BEDTOOLS_MERGE.out.bed        // channel: [ val(meta), [ bed ] ]
    filtered_bed      = AWK.out.file                  // channel: [ val(meta), [ bed ] ]
    versions          = ch_versions                   // channel: [ versions.yml       ]
}
