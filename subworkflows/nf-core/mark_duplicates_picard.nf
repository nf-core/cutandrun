/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'

workflow MARK_DUPLICATES_PICARD {
    take:
    bam            // channel: [ val(meta), [ bam ] ]
    process_target //boolean

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
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
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
        ch_bam      = PICARD_MARKDUPLICATES.out.bam.mix ( ch_split.target )
        metrics     = PICARD_MARKDUPLICATES.out.metrics
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    }
    //ch_bam | view

    /*
    * Index BAM file
    */
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Join bam/bai
    ch_bam
        .map { row -> [row[0].id, row ].flatten()}
        .join ( SAMTOOLS_INDEX.out.bai.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
    //ch_bam_bai | view

    /*
    * Run samtools stats, flagstat and idxstats
    */
    BAM_STATS_SAMTOOLS ( ch_bam_bai )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = ch_bam                            // channel: [ val(meta), [ bam ] ]
    metrics                                      // channel: [ val(meta), [ metrics ] ]

    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
