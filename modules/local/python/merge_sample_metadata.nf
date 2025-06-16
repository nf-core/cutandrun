process MERGE_SAMPLE_METADATA {
    label 'process_single'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path metadata

    output:
    path '*.csv'             , emit: csv
    path  "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ext             = task.ext.output ?: 'csv'
    def output          = task.ext.output ?: 'NO_NAME.csv'
    def id_parse_string = task.ext.id_parse_string ?: ''

    """
    reports.py merge_samples \\
        --metadata "*.${ext}" \\
        --id_parse_string $id_parse_string \\
        --output $output \\
        --log log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
