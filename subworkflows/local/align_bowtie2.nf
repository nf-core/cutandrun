/*
 * Alignment with BOWTIE2
 */

include { BOWTIE2_ALIGN as BOWTIE2_TARGET_ALIGN                      } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN as BOWTIE2_SPIKEIN_ALIGN                     } from '../../modules/nf-core/bowtie2/align/main'
include { BAM_SORT_STATS_SAMTOOLS                                    } from '../nf-core/bam_sort_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_SPIKEIN } from '../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_BOWTIE2 {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    index         // channel: [ val(meta), [ index ] ]
    spikein_index // channel: [ val(meta), [ index ] ]
    fasta         // channel: [ val(meta), fasta ]
    spikein_fasta // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    /*
     * Map reads with BOWTIE2 to target genome
     */
    ch_index = index.map { [[id:it.baseName], it] }
    BOWTIE2_TARGET_ALIGN (
        reads,
        ch_index.collect{ it[1] },
        params.save_unaligned,
        false
    )
    ch_versions = ch_versions.mix(BOWTIE2_TARGET_ALIGN.out.versions)

    /*
     * Map reads with BOWTIE2 to spike-in genome
     */
    ch_spikein_index = spikein_index.map { [[id:it.baseName], it] }
    BOWTIE2_SPIKEIN_ALIGN (
        reads,
        ch_spikein_index.collect{ it[1] },
        params.save_unaligned,
        false
    )
    ch_versions = ch_versions.mix(BOWTIE2_SPIKEIN_ALIGN.out.versions)

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_STATS_SAMTOOLS ( BOWTIE2_TARGET_ALIGN.out.aligned, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    BAM_SORT_STATS_SAMTOOLS_SPIKEIN ( BOWTIE2_SPIKEIN_ALIGN.out.aligned, spikein_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.versions)

    emit:
    versions             = ch_versions                                  // channel: [ versions.yml ]

    orig_bam             = BOWTIE2_TARGET_ALIGN.out.aligned             // channel: [ val(meta), bam ]
    orig_spikein_bam     = BOWTIE2_SPIKEIN_ALIGN.out.aligned            // channel: [ val(meta), bam ]

    bowtie2_log          = BOWTIE2_TARGET_ALIGN.out.log                 // channel: [ val(meta), log_final ]
    bowtie2_spikein_log  = BOWTIE2_SPIKEIN_ALIGN.out.log                // channel: [ val(meta), log_final ]

    bam                  = BAM_SORT_STATS_SAMTOOLS.out.bam              // channel: [ val(meta), [ bam ] ]
    bai                  = BAM_SORT_STATS_SAMTOOLS.out.bai              // channel: [ val(meta), [ bai ] ]
    stats                = BAM_SORT_STATS_SAMTOOLS.out.stats            // channel: [ val(meta), [ stats ] ]
    flagstat             = BAM_SORT_STATS_SAMTOOLS.out.flagstat         // channel: [ val(meta), [ flagstat ] ]
    idxstats             = BAM_SORT_STATS_SAMTOOLS.out.idxstats         // channel: [ val(meta), [ idxstats ] ]

    spikein_bam          = BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.bam      // channel: [ val(meta), [ bam ] ]
    spikein_bai          = BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.bai      // channel: [ val(meta), [ bai ] ]
    spikein_stats        = BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.stats    // channel: [ val(meta), [ stats ] ]
    spikein_flagstat     = BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    spikein_idxstats     = BAM_SORT_STATS_SAMTOOLS_SPIKEIN.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
