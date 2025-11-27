/*
 * Extract fragments from a BAM file into a bedfile format
*/

include { AWK               } from '../../modules/local/linux/awk'
include { BEDTOOLS_BAMTOBED } from "../../modules/nf-core/bedtools/bamtobed"
include { CUT               } from '../../modules/local/linux/cut'
include { SAMTOOLS_SORT     } from "../../modules/nf-core/samtools/sort"

workflow EXTRACT_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Sort reads by name for bamtobed
    */
    SAMTOOLS_SORT(bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
    * MODULE: Convert BAM file to paired-end bed format
    */
    BEDTOOLS_BAMTOBED(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)
    // BEDTOOLS_BAMTOBED.out.bed | view

    /*
    * MODULE: Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    */
    AWK(BEDTOOLS_BAMTOBED.out.bed)
    ch_versions = ch_versions.mix(AWK.out.versions)

    /*
    * MODULE: Only extract the fragment related columns
    */
    CUT(AWK.out.file)
    ch_versions = ch_versions.mix(CUT.out.versions)

    emit:
    bed      = CUT.out.file // channel: [ val(meta), [ bed ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
