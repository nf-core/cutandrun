process FILTER_BAMS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam)         , optional: true, emit: bam
    path  "versions.yml"               , emit: versions


    when:
    num_alignments > 1000

    script:
    """
    num_alignments=\$(samtools view -c ${bam})
    echo "Number of alignments in ${bam}: \${num_alignments}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
