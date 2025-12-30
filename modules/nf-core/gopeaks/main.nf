process GOPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gopeaks=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gopeaks:1.0.0--h9ee0642_0' :
        'biocontainers/gopeaks:1.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(control_bam)

    output:
    tuple val(meta), path("*_peaks.bed"), emit: peaks
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam ? "-c $control_bam" : ''
    """
    gopeaks \\
        -b $bam \\
        $control \\
        -o ${prefix}_peaks.bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gopeaks: \$(gopeaks --version 2>&1 | sed 's/gopeaks version //')
    END_VERSIONS
    """
}
