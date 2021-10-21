/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

params.markduplicates_options = [:]
params.samtools_options       = [:]
params.control_only           = false

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main' addParams( options: params.markduplicates_options )
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'        addParams( options: params.samtools_options       )
include { BAM_STATS_SAMTOOLS    } from './bam_stats_samtools'                                     addParams( options: params.samtools_options       )

workflow MARK_DUPLICATES_PICARD {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    /*
    * Picard MarkDuplicates
    */
    out_bam = Channel.empty()
    metrics = Channel.empty()
    version = Channel.empty()
    if( !params.control_only ) {
        PICARD_MARKDUPLICATES ( bam )
        out_bam = PICARD_MARKDUPLICATES.out.bam
        metrics = PICARD_MARKDUPLICATES.out.metrics
        version = PICARD_MARKDUPLICATES.out.version
    }
    else { // Split out non igg files and run only on these
        bam.branch { it ->
            target: it[0].group != 'igg'
            control: it[0].group == 'igg'
        }
        .set { ch_split }

        PICARD_MARKDUPLICATES ( ch_split.control )
        out_bam = PICARD_MARKDUPLICATES.out.bam
        metrics = PICARD_MARKDUPLICATES.out.metrics
        version = PICARD_MARKDUPLICATES.out.version

        out_bam = out_bam.mix ( ch_split.target )
    }
    //out_bam | view

    /*
    * Index BAM file
    */
    SAMTOOLS_INDEX ( out_bam )

    // Join bam/bai
    out_bam
        .map { row -> [row[0].id, row ].flatten()}
        .join ( SAMTOOLS_INDEX.out.bai.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
    //ch_bam_bai | view

    /*
    * Run samtools stats, flagstat and idxstats
    */
    BAM_STATS_SAMTOOLS ( ch_bam_bai )

    emit:
    bam              = out_bam                           // channel: [ val(meta), [ bam ] ]
    metrics                                              // channel: [ val(meta), [ metrics ] ]
    picard_version   = version                           // path: *.version.txt

    bai              = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), [ bai ] ]
    stats            = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    samtools_version = SAMTOOLS_INDEX.out.version        // path: *.version.txt
}
