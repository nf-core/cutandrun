/*
 * Calculate bed fragments from bam file
 */

include { SAMTOOLS_VIEW      } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'
include { BEDTOOLS_BAMTOBED  } from '../../modules/nf-core/modules/bedtools/bamtobed/main'
include { AWK                } from '../../modules/local/linux/awk'
include { CUT                } from '../../modules/local/linux/cut'

workflow CALCULATE_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    ch_versions = Channel.empty()

    /*
     * Filter BAM file
     */
    SAMTOOLS_VIEW ( bam, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    /*
    * Sort BAM file
    */
    SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Convert to bed file
    BEDTOOLS_BAMTOBED ( SAMTOOLS_SORT.out.bam  )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    // Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    AWK ( BEDTOOLS_BAMTOBED.out.bed  )
    ch_versions = ch_versions.mix(AWK.out.versions)

    // Only extract the fragment related columns
    CUT ( AWK.out.file )

    emit:
    bed              = CUT.out.file                   // channel: [ val(meta), [ bed ] ]
    bam              = SAMTOOLS_SORT.out.bam          // channel: [ val(meta), [ bam ] ]
    versions         = ch_versions                    // channel: [ versions.yml ]
}
