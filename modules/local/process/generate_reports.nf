include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GENERATE_REPORTS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    //container "quay.io/biocontainers/python:3.8.3"
    //container 'quay.io/biocontainers/pybda:0.1.0--pyh5ca1d4c_0'
    container "luslab/cutandrun-dev-reporting:latest"

    input:
    path meta_data
    path raw_fragments
    
    output:
    path '*.pdf', emit: pdf
    path '*.csv', emit: csv
    path '*.png', emit: png

    script:  // This script is bundled with the pipeline, in nf-core/cutandrun/bin/
    """
    $baseDir/bin/python/reporting/main.py genimg --meta $meta_data --raw_frag "*_raw.csv" --output . --log log.txt
    """
}