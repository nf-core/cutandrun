/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main'
include { BAM_SORT_SAMTOOLS     } from './bam_sort_samtools'

workflow MARK_DUPLICATES_PICARD {
    take:
    bam            // channel: [ val(meta), [ bam ] ]
    process_target // boolean

    main:
    /*
    * Picard MarkDuplicates
    */
    ch_bam      = Channel.empty()
    metrics     = Channel.empty()
    ch_versions = Channel.empty()
    if( process_target ) {
        PICARD_MARKDUPLICATES ( bam )
        ch_bam      = PICARD_MARKDUPLICATES.out.bam
        metrics     = PICARD_MARKDUPLICATES.out.metrics
        ch_versions = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions )
    }
    else { // Split out control files and run only on these
        bam.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_split }
        //ch_split.target | view
        //ch_split.control | view

        PICARD_MARKDUPLICATES ( ch_split.control )

        // Prevents issues with resume with the branch elements coming in the wrong order
        ch_sorted_targets = ch_split.target
            .toSortedList { row -> row[0].id }
            .flatMap()

        ch_sorted_controls = PICARD_MARKDUPLICATES.out.bam
            .toSortedList { row -> row[0].id }
            .flatMap()

        ch_bam      = ch_sorted_targets.concat ( ch_sorted_controls )
        metrics     = PICARD_MARKDUPLICATES.out.metrics
        ch_versions = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions )
    }
    //ch_bam | view

    /*
    * WORKFLOW: Re sort and index all the bam files + calculate stats
    */
    BAM_SORT_SAMTOOLS ( 
        ch_bam 
    )

    emit:
    bam      = BAM_SORT_SAMTOOLS.out.bam        // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_SAMTOOLS.out.bai        // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    metrics                                     // channel: [ val(meta), [ metrics ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}
