process CUT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.cut.*"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext     = task.ext.ext ?: 'txt'
    def command = task.ext.command ?: ''

    """
    cut $args $input $command > ${prefix}.cut.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version | head -n 1 | awk '{print \$4;}')
    END_VERSIONS
    """
}
