process CREATE_MQC {
    tag "$meta.id"
    label 'process_ultralow'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(files)
    path header

    output:
    path("*_mqc*"), emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    def ext    = task.ext.ext ?: "tsv" 
    """
    echo "total_peaks" > id.txt
    paste id.txt ${files.join(' ')} > data.txt
    cat $header data.txt > ${prefix}_mqc.${ext}
    """
}
