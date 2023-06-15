process FRAG_LEN_HIST {
    label 'process_medium'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path raw_fragments
    path frag_len_header_multiqc

    output:
    path '*frag_len_mqc.yml', emit: frag_len_mqc
    path  "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    calc_frag_hist.py \\
        --frag_path "*len.txt" \\
        --output frag_len_hist.txt

    if [ -f "frag_len_hist.txt" ]; then
        cat $frag_len_header_multiqc frag_len_hist.txt > frag_len_mqc.yml
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        numpy: \$(python -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        seaborn: \$(python -c 'import seaborn; print(seaborn.__version__)')
    END_VERSIONS
    """
}
