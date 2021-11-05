/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

params.options = [:]
params.samtools_sort_options = [:]

include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'  addParams( options: params.samtools_sort_options )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: params.options               )
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools'                              addParams( options: params.options               )

workflow BAM_SORT_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
    /*
    * SORT BAM file
    */
    SAMTOOLS_SORT      ( ch_bam )

    /*
    * Index BAM file
    */
    SAMTOOLS_INDEX     ( SAMTOOLS_SORT.out.bam )

    // Join bam/bai
    ch_bam_sample_id = SAMTOOLS_SORT.out.bam.map  { row -> [row[0].id, row] }
    ch_bai_sample_id = SAMTOOLS_INDEX.out.bai.map { row -> [row[0].id, row] }
    ch_bam_bai = ch_bam_sample_id.join(ch_bai_sample_id, by: [0]).map {row -> [row[1][0], row[1][1], row[2][1]]}

    /*
    * Run samtools stats, flagstat and idxstats
    */
    BAM_STATS_SAMTOOLS ( ch_bam_bai )

    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    version  = SAMTOOLS_SORT.out.version       // path: *.version.txt
}
