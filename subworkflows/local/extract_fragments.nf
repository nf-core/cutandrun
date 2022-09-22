/*
 * Extract fragments from a BAM file into a bedfile format
*/

include { BEDTOOLS_BAMTOBED  } from "../../modules/nf-core/modules/bedtools/bamtobed/main.nf"
include { AWK                } from '../../modules/local/linux/awk'
include { CUT                } from '../../modules/local/linux/cut'

workflow EXTRACT_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Convert BAM file to paired-end bed format
    */
    BEDTOOLS_BAMTOBED(
        bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)
    // BEDTOOLS_BAMTOBED.out.bed | view

    /*
    * MODULE: Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    */
    AWK (
        BEDTOOLS_BAMTOBED.out.bed  
    )
    ch_versions = ch_versions.mix(AWK.out.versions)

    /*
    * MODULE: Only extract the fragment related columns
    */
    CUT (
        AWK.out.file
    )

    emit:
    bed      = CUT.out.file // channel: [ val(meta), [ bed ] ]
    versions = ch_versions  // channel: [ versions.yml ]
}
