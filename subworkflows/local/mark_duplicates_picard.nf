/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { PICARD_MARKDUPLICATES   } from '../../modules/nf-core/picard/markduplicates/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow MARK_DUPLICATES_PICARD {
    take:
    bam            // channel: [ val(meta), [ bam ] ]
    bai            // channel: [ val(meta), [ bai ] ]
    process_target // boolean
    fasta          // channel: [ val(meta), fasta ]
    fai            // channel: [ val(meta), fai ]

    main:
    // Init
    ch_bam      = Channel.empty()
    ch_bai      = Channel.empty()
    ch_metrics  = Channel.empty()
    ch_stats    = Channel.empty()
    ch_flagstat = Channel.empty()
    ch_idxstats = Channel.empty()
    ch_versions = Channel.empty()

    if( process_target ) {
        PICARD_MARKDUPLICATES (
            bam,
            fasta,
            fai
        )
        ch_metrics  = PICARD_MARKDUPLICATES.out.metrics
        ch_versions = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions )

        /*
        * WORKFLOW: Re sort and index all the bam files + calculate stats
        */
        BAM_SORT_STATS_SAMTOOLS (
            PICARD_MARKDUPLICATES.out.bam,
            fasta
        )
        ch_bam      = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_bai      = BAM_SORT_STATS_SAMTOOLS.out.bai
        ch_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
        ch_versions = ch_versions.mix( BAM_SORT_STATS_SAMTOOLS.out.versions )
    }
    else { // Split out control files and run only on these
        bam.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_split_bam }
        //ch_split.target | view
        //ch_split.control | view

        bai.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_split_bai }

        PICARD_MARKDUPLICATES (
            ch_split_bam.control,
            fasta,
            fai
        )
        ch_metrics  = PICARD_MARKDUPLICATES.out.metrics
        ch_versions = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions )

        /*
        * WORKFLOW: Re sort and index all the bam files + calculate stats
        */
        BAM_SORT_STATS_SAMTOOLS (
            PICARD_MARKDUPLICATES.out.bam,
            fasta
        )
        ch_bam      = BAM_SORT_STATS_SAMTOOLS.out.bam.mix(ch_split_bam.target)
        ch_bai      = BAM_SORT_STATS_SAMTOOLS.out.bai.mix(ch_split_bai.target)
        ch_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
        ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
        ch_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
        ch_versions = ch_versions.mix( BAM_SORT_STATS_SAMTOOLS.out.versions )
    }

    emit:
    bam      = ch_bam      // channel: [ val(meta), [ bam ] ]
    bai      = ch_bai      // channel: [ val(meta), [ bai ] ]
    stats    = ch_stats    // channel: [ val(meta), [ stats ] ]
    flagstat = ch_flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = ch_idxstats // channel: [ val(meta), [ idxstats ] ]
    metrics  = ch_metrics  // channel: [ val(meta), [ metrics ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
