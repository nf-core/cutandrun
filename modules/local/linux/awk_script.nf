process AWK_SCRIPT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'quay.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(input)
    path script

    output:
    tuple val(meta), path("*.awk.txt"), emit: file
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    awk $args -f $script $input > ${prefix}.awk.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version | head -n 1 | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
