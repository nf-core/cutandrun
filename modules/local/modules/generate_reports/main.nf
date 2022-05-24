process GENERATE_REPORTS {
    label 'process_ultralow'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas=1.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path meta_data
    path meta_data_ctrl
    path raw_fragments
    path bed_fragments
    path seacr_beds
    path frag_len_header_multiqc

    output:
    path '*.pdf'             , emit: pdf
    path '*.csv'             , emit: csv
    path '*.png'             , emit: png
    path '*frag_len_mqc.yml', emit: frag_len_multiqc
    path  "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def meta_data_resolved = meta_data ? meta_data : meta_data_ctrl

    """
    reporting.py gen_reports \\
        --meta $meta_data_resolved \\
        --meta_ctrl $meta_data_ctrl \\
        --raw_frag "*.frags.len.txt" \\
        --bin_frag "*bin500.awk.bed" \\
        --seacr_bed "*bed*.bed" \\
        --output . \\
        --log log.txt

    if [ -f "03_03_frag_len_mqc.txt" ]; then
        cat $frag_len_header_multiqc 03_03_frag_len_mqc.txt > frag_len_mqc.yml
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """

}
