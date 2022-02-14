/*
 * Annotate the pipeline meta data with the columns from a csv file
 * generated from processing a report text file with an awk script
 */

include { AWK_SCRIPT } from '../../modules/local/linux/awk_script'
include { AWK }        from '../../modules/local/linux/awk'

workflow ANNOTATE_META_AWK {
    take:
    passthrough
    report
    script
    meta_suffix // string
    meta_prefix // string
    script_mode // bool

    main:
    ch_versions = Channel.empty()

    // Strip out the sample id from the meta in the passthrough
    ch_paths = passthrough.map { row -> [row[0].id, row[0], row[1..-1]].flatten() }

    ch_annotated_meta = Channel.empty()
    // Can run awk in script mode with a file from assets or with a setup of command line args
    if(script_mode) {
        AWK_SCRIPT ( report, script )
        ch_versions = ch_versions.mix(AWK_SCRIPT.out.versions)

        AWK_SCRIPT.out.file
            .splitCsv(header:true)
            .map { row ->
                new_meta = [:]
                row[1].each{ k, v -> new_meta.put(meta_prefix + k + meta_suffix, v) }
                [row[0].id, new_meta]
            }
            .join ( ch_paths )
            .map { row -> [ row[2] << row[1], row[3..-1] ] }
            .set { ch_annotated_meta }
    }
    else {
        AWK ( report )
        ch_versions = ch_versions.mix(AWK.out.versions)

        AWK.out.file
            .splitCsv(header:true)
            .map { row ->
                new_meta = [:]
                row[1].each{ k, v -> new_meta.put(meta_prefix + k + meta_suffix, v) }
                [row[0].id, new_meta]
            }
            .join ( ch_paths )
            .map { row -> [ row[2] << row[1], row[3..-1] ] }
            .set { ch_annotated_meta }
    }

    emit:
    output     = ch_annotated_meta // channel: [ val(annotated_meta), [ passthrough ] ]
    versions   = ch_versions       // channel: [ versions.yml ]
}
