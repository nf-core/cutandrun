/*
 * Annotate the pipeline meta data with a csv file 
 */

params.options     = [:]
params.meta_suffix = ""
params.meta_prefix = ""

workflow ANNOTATE_META_CSV {
    take: passthrough
    take: reports
    main:

    main:
    // Strip out the sample id from the meta in the passthrough
    ch_paths = passthrough.map { row -> [row[0].id, row[0], row[1..-1]].flatten() }

    reports.splitCsv(header:true)
       .map { row ->
            new_meta = [:]
            row[1].each{ k, v -> new_meta.put(params.meta_prefix + k + params.meta_suffix, v) }
            [row[0].id, new_meta]
        }
        .join ( ch_paths )
        .map { row -> [ row[2] << row[1], row[3..-1] ] } 
        .set { ch_annotated_meta }

    emit:
    output = ch_annotated_meta // channel: [ val(annotated_meta), [ passthrough ] ]
}
