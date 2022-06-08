/*
 * Annotate the pipeline meta data with a csv file
 */

workflow ANNOTATE_META_CSV {
    take:
    passthrough // channel
    reports     // file
    meta_suffix // string
    meta_prefix // string

    main:

    main:
    // Strip out the sample id from the meta in the passthrough
    ch_paths = passthrough.map { row -> [row[0].id, row[0], row[1..-1]].flatten() }

    reports.splitCsv(header:true)
        .map { row ->
            new_meta = [:]
            row[1].each{ k, v -> new_meta.put(meta_prefix + k + meta_suffix, v) }
            [row[0].id, new_meta]
        }
        .join ( ch_paths )
        .map { row -> [ row[2] << row[1], row[3..-1] ] }
        .set { ch_annotated_meta }

    emit:
    output = ch_annotated_meta // channel: [ val(annotated_meta), [ passthrough ] ]
}
