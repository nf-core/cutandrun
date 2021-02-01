include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GENERATE_REPORTS {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    //container "quay.io/biocontainers/python:3.8.3"
    container 'quay.io/biocontainers/pybda:0.1.0--pyh5ca1d4c_0'
    
    output:
    path '*.pdf'


    script:  // This script is bundled with the pipeline, in nf-core/cutandrun/bin/
    """
    gen_pdf.py
    """
}