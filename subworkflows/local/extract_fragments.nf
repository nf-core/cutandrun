/*
 * Calculate bed fragments from bam file
 */

include { BEDTOOLS_BAMTOBED   } from '../../modules/nf-core/modules/bedtools/bamtobed/main'
include { AWK                 } from '../../modules/local/linux/awk'
include { CUT                 } from '../../modules/local/linux/cut'
include { AWK as AWK_FRAG_BIN } from "../../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW } from "../../modules/local/samtools_custom_view"

workflow EXTRACT_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    bai // channel: [ val(meta), [ bai ] ]

    main:
    ch_versions = Channel.empty()

    /*
    * CHANNEL: Filter bams for target only
    */
    bam.filter { it -> it[0].is_control == false }
    .set { ch_bam_target }
    //ch_bam_target | view

    /*
    * CHANNEL: Filter bais for target only
    */
    bai.filter { it -> it[0].is_control == false }
    .set { ch_bai_target }
    //ch_bai_target | view

    /*
    * MODULE: Convert to bed file
    */
    BEDTOOLS_BAMTOBED ( ch_bam_target  )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    /*
    * MODULE: Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    */
    AWK ( BEDTOOLS_BAMTOBED.out.bed  )
    ch_versions = ch_versions.mix(AWK.out.versions)

    /*
    * MODULE: Only extract the fragment related columns
    */
    CUT ( AWK.out.file )
    ch_versions = ch_versions.mix(CUT.out.versions)

    /*
    * MODULE: Bin the fragments into 500bp bins ready for downstream reporting
    */
    AWK_FRAG_BIN(
        CUT.out.file
    )
    ch_versions = ch_versions.mix(AWK_FRAG_BIN.out.versions)
    //AWK_FRAG_BIN.out.file | view

    /*
    * CHANNEL: Combine bam and bai files on id
    */
    bam.map { row -> [row[0].id, row ].flatten()}
    .join ( ch_bai_target.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [row[1], row[2], row[4]] }
    .set { ch_bam_bai }
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bam_bai | view

    /*
    * MODULE: Calculate fragment lengths
    */
    SAMTOOLS_CUSTOMVIEW (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CUSTOMVIEW.out.versions)
    //SAMTOOLS_CUSTOMVIEW.out.tsv | view

    emit:
    bed          = CUT.out.file                // channel: [ val(meta), [ bed ] ]
    binned_frags = AWK_FRAG_BIN.out.file       // channel: [ val(meta), [ bed ] ]
    frag_lengths = SAMTOOLS_CUSTOMVIEW.out.tsv // channel: [ val(meta), [ tsv ] ]
    versions     = ch_versions                 // channel: [ versions.yml ]
}
