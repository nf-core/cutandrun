process EPIC2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::epic2=0.0.52"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epic2:0.0.52--py39h14c64f4_0' :
        'biocontainers/epic2:0.0.52--py39h14c64f4_0' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    val genome

    output:
    tuple val(meta), path("*.peaks"), emit: peaks
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epic2 \\
        --treatment $treatment_bam \\
        --control $control_bam \\
        --genome $genome \\
        --output ${prefix}.peaks \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epic2: \$(epic2 --version 2>&1 | sed 's/epic2 //')
    END_VERSIONS
    """
}
