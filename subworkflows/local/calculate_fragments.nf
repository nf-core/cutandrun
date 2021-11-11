/*
 * Calculate bed fragments from bam file
 */

params.samtools_options      = [:]
params.samtools_view_options = [:]
params.samtools_sort_options = [:]
params.bamtobed_options      = [:]
params.awk_options           = [:]
params.cut_options           = [:]

include { SAMTOOLS_VIEW      } from "../../modules/nf-core/modules/samtools/view/main"     addParams( options: params.samtools_view_options )
include { SAMTOOLS_SORT      } from "../../modules/nf-core/modules/samtools/sort/main"     addParams( options: params.samtools_sort_options )
include { BEDTOOLS_BAMTOBED  } from "../../modules/nf-core/modules/bedtools/bamtobed/main" addParams( options: params.bamtobed_options      )
include { AWK                } from "../../modules/local/linux/awk"                        addParams( options: params.awk_options           )
include { CUT                } from "../../modules/local/linux/cut"                        addParams( options: params.cut_options           )

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
