/*
 * Annotate the pipeline meta data with the columns from a csv file
   generated from processing a report text file with an awk script
 */

params.options = [:]

include { AWK } from '../process/awk' addParams( options: params.options )

workflow ANNOTATE_META {
    take: input
    take: report
    take: script
    main:
    
    main:
    // Strip sample id and paths only from input
    ch_paths = input.map { row -> [row[0].id, row[1..-1]].flatten() }

    // Create csv file from the report channel using the awk script
    AWK ( report, script )

    // Annotate the input meta data with the csv values generated from 
    // the awk script
    AWK.out.file
        .splitCsv(header:true)
        .map { row -> [ row[0].id, row[0] << row[1] ] }
        .join ( ch_paths )
        .map { row -> row[1..-1] }
        .set { ch_annotated_meta }

    emit:
    output = ch_annotated_meta // channel: [ val(meta), [ input ] ]
}
