process AWK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.awk.*"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext      = task.ext.ext ?: 'txt'
    def command  = task.ext.command ?: ''
    def command2 = task.ext.command2 ?: ''

    """
    awk $args $command $input $command2 > ${prefix}.awk.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>/dev/null | head -n 1 | awk '{split(\$0,a,","); print a[1];}' | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
