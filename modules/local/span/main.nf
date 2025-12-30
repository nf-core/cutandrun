process SPAN {
    tag "$meta.id"
    label 'process_high_memory'

    conda "bioconda::span=0.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/span:0.14.1--hdfd78af_0' :
        'biocontainers/span:0.14.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    path chrom_sizes

    output:
    tuple val(meta), path("*.peak"), emit: peaks
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = control_bam.name != 'NO_CONTROL' ? "-c $control_bam" : ''
    def memory = task.memory ? task.memory.toGiga() : 8
    """
    span \\
        analyze \\
        -t $treatment_bam \\
        $control \\
        --chrom.sizes $chrom_sizes \\
        -o ${prefix}.peak \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        span: \$(span --version 2>&1 | head -n1 | sed 's/SPAN //g')
    END_VERSIONS
    """
}
