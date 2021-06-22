include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options      = initOptions(params.options)

process PLOT_CONSENSUS_PEAKS {
    // tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }


    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::numpy=1.20.* conda-forge::pandas=1.2.* conda-forge::upsetplot=0.4.4" : null)
    container "luslab/cutandrun-dev-plot-consensus-peaks:latest"

    input:
    // tuple val(meta), path(consensus_peaks)
    path(consensus_peaks)

    output:
    path ("*.pdf"), optional:true, emit: pdf
    path '*.version.txt', emit: version

    script:
    """
    consensus_peaks.py \\
        --peaks $consensus_peaks \\
        --outpath .

    python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\" > python.version.txt
    """

}
