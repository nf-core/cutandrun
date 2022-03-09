process PLOT_CONSENSUS_PEAKS {
    label 'process_min'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::numpy=1.20.* conda-forge::pandas=1.2.* conda-forge::upsetplot=0.4.4" : null)
    container "luslab/cutandrun-dev-plot-consensus-peaks:latest"

    input:
    path(consensus_peaks)

    output:
    path ("*.pdf")      , optional:true, emit: pdf
    path  "versions.yml", emit: versions

    script:
    """
    consensus_peaks.py \\
        --peaks "*.peaks.bed" \\
        --outpath .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """

}
