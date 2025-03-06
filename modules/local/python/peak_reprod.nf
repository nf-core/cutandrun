process CALCULATE_PEAK_REPROD {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.8.3 conda-forge::dask=2021.9.1 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(bed)
    path peak_reprod_header_multiqc

    output:
    tuple val(meta), path("*peak_repro.tsv"), emit: tsv
    path "*_mqc.tsv"                        , emit: mqc
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    peak_reproducibility.py \\
        --sample_id $meta.id \\
        --intersect $bed \\
        --threads ${task.cpus} \\
        --outpath .

    cat $peak_reprod_header_multiqc *peak_repro.tsv > ${prefix}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        dask: \$(python -c 'import dask; print(dask.__version__)')
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
