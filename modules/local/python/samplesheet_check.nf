// Import generic module functions
include { initOptions; saveFiles } from '../common/functions'

params.options = [:]
options        = initOptions(params.options)

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet

    output:
    path '*.csv'


    script:  // This script is bundled with the pipeline, in nf-core/cutandrun/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv $params.igg_control
    """
}
