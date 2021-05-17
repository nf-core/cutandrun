/*
 * Run bam files through samtools view and reindex and calc stats
 */

params.samtools_view_options = [:]
params.samtools_options      = [:]

include { SAMTOOLS_VIEW      } from '../../modules/nf-core/software/samtools/view/main'  addParams( options: params.samtools_view_options )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/software/samtools/index/main' addParams( options: params.samtools_options      )
include { BAM_STATS_SAMTOOLS } from '../nf-core/bam_stats_samtools'                      addParams( options: params.samtools_options      )

workflow SAMTOOLS_VIEW_SORT_STATS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    /*
     * Filter using samtools view
     */
    SAMTOOLS_VIEW ( bam )

    /*
     * Index BAM file and run samtools stats, flagstat and idxstats
     */
    SAMTOOLS_INDEX     ( SAMTOOLS_VIEW.out.bam )
    BAM_STATS_SAMTOOLS ( SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0]) )

    emit:
    bam              = SAMTOOLS_VIEW.out.bam           // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_INDEX.out.version      //    path: *.version.txt
}
