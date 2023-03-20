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
    tuple val(meta), path(bam), stdout  , emit: bam

    when:
    task.ext.when == null || task.ext.when 

    script:
    """
    if [[ \$(samtools view -c $bam) -ge 1000 ]];
    then
        echo "1"
    else
        echo "0"
    fi
    """
}
