process DEEPTOOLS_PLOT_FINGERPRINT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*.raw.txt'                  , emit: fingerprint
    path '*.{txt,pdf}'                , emit: quality_metrics
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotFingerprint $args \\
        --bamfiles ${bam} \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --labels $prefix \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --skipZeros \\
        --numberOfProcessors $task.cpus
        
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotFingerprint --version | sed -e "s/plotFingerprint //g")
    END_VERSIONS
    """
}