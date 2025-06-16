/*
 * Read QC, read trimming and post trim QC
 */

include { FASTQC     } from '../../modules/nf-core/fastqc/main'
include { TRIMGALORE } from '../../modules/local/for_patch/trimgalore/main'

workflow FASTQC_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false

    main:
    ch_versions = Channel.empty()

    fastqc_html     = Channel.empty()
    fastqc_zip      = Channel.empty()
    if (!skip_fastqc) {
        FASTQC ( reads ).html.set { fastqc_html }
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions)
    }

    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    ch_output_reads = reads
    if (!skip_trimming) {
        TRIMGALORE (
            reads
        )
        ch_output_reads = TRIMGALORE.out.reads
        trim_html       = TRIMGALORE.out.html
        trim_zip        = TRIMGALORE.out.zip
        trim_log        = TRIMGALORE.out.log
        ch_versions     = ch_versions.mix(TRIMGALORE.out.versions)
    }

    emit:
    reads    = ch_output_reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions     // channel: [ versions.yml ]

    fastqc_html                // channel: [ val(meta), [ html ] ]
    fastqc_zip                 // channel: [ val(meta), [ zip ] ]

    trim_html                  // channel: [ val(meta), [ html ] ]
    trim_zip                   // channel: [ val(meta), [ zip ] ]
    trim_log                   // channel: [ val(meta), [ txt ] ]
}
