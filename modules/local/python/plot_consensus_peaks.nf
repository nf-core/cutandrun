process PLOT_CONSENSUS_PEAKS {
    label 'process_single'

    conda "conda-forge::python=3.8.3 conda-forge::numpy=1.20.* conda-forge::pandas=1.2.* conda-forge::upsetplot=0.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path(consensus_peaks)

    output:
    path ("*.pdf")      , optional:true, emit: pdf
    path  "versions.yml", emit: versions

    script:
    """
    plot_consensus_peaks.py \\
        --peaks "*.peaks.bed" \\
        --outpath .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        upsetplot: \$(python -c 'import upsetplot; print(upsetplot.__version__)')
    END_VERSIONS
    """
}
