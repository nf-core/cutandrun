process CUT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'quay.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.cut.*"), emit: file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext     = task.ext.ext ?: 'txt'
    def command = task.ext.command ?: ''

    """
    cut $args $input $command > ${prefix}.cut.${ext}
    """
}
