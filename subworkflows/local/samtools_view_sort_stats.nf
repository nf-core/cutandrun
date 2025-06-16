/*
 * Run bam files through samtools view and reindex and calc stats
 */

include { SAMTOOLS_VIEW      } from '../../modules/local/for_patch/samtools/view/main'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../nf-core/bam_stats_samtools/main'

workflow SAMTOOLS_VIEW_SORT_STATS {
    take:
    bam     // channel: [ val(meta), [ bam ] ]
    regions // channel: [ regions ]
    fasta   // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()
    /*
     * Filter BAM file
     */
    SAMTOOLS_VIEW ( bam.map{ row -> [ row[0], row[1], [], ] }, [], regions )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    /*
    * Sort BAM file
    */
    SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    /*
     * Index BAM file and run samtools stats, flagstat and idxstats
     */
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Join bam/bai
    ch_bam_sample_id = SAMTOOLS_SORT.out.bam.map  { row -> [row[0].id, row] }
    ch_bai_sample_id = SAMTOOLS_INDEX.out.bai.map { row -> [row[0].id, row] }
    ch_bam_bai = ch_bam_sample_id.join(ch_bai_sample_id, by: [0]).map {row -> [row[1][0], row[1][1], row[2][1]]}

    BAM_STATS_SAMTOOLS ( ch_bam_bai, fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions.first())

    emit:
    bam      = SAMTOOLS_SORT.out.bam            // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai           // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats     // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats  // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
