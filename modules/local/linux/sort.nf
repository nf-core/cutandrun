process SORT {
    tag "$meta.id"
    label 'process_ultralow'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.sort.*"), emit: file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext     = task.ext.ext ?: 'txt'

    // Option for multiple files to sort
    String input_files = input.join(" ")

    """
    sort -T '.' $args $input_files > ${prefix}.sort.${ext}
    """
}
