/*
 * Annotate the pipeline meta data with the columns from a csv file
   generated from processing a report text file with an awk script
   
   TODO: input is a passthrough and that is the thing that is annotated
    */

params.options = [:]
params.meta_suffix = ''
params.meta_prefix = ''
params.script_mode = false

include { AWK_SCRIPT } from '../process/awk_script' addParams( options: params.options )
include { AWK } from '../process/awk' addParams( options: params.options )

workflow ANNOTATE_META_AWK {
    take: passthrough
    take: report
    take: script
    main:
    
    main:
    // Strip out the sample id from the meta in the passthrough
    ch_paths = passthrough.map { row -> [row[0].id, row[0], row[1..-1]].flatten() }

    ch_annotated_meta = Channel.empty()
    if(params.script_mode) {
        AWK_SCRIPT ( report, script )

        AWK_SCRIPT.out.file
            .splitCsv(header:true)
            .map { row -> 
                new_meta = [:]
                row[1].each{ k, v -> new_meta.put(params.meta_prefix + k + params.meta_suffix, v) }
                [row[0].id, new_meta]
            }
            .join ( ch_paths )
            .map { row -> [ row[2] << row[1], row[3..-1] ] }
            .set { ch_annotated_meta }
    }
    else {
        AWK ( report )

        AWK.out.file
            .splitCsv(header:true)
            .map { row -> 
                new_meta = [:]
                row[1].each{ k, v -> new_meta.put(params.meta_prefix + k + params.meta_suffix, v) }
                [row[0].id, new_meta]
            }
            .join ( ch_paths )
            .map { row -> [ row[2] << row[1], row[3..-1] ] }
            .set { ch_annotated_meta }
    }

    emit:
    output = ch_annotated_meta // channel: [ val(annotated_meta), [ passthrough ] ]
}
