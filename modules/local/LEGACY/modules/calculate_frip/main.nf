process CALCULATE_FRIP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::deeptools=3.5.* bioconda::pysam=0.17.*" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path('frips.csv'), emit: frips
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ext = bed.getExtension()
    """
    frip.py \\
        --bams "*.bam" \\
        --peaks "*.${ext}" \\
        --threads ${task.cpus} \\
        --outpath .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        deeptools: \$(deeptools --version | sed -e "s/deeptools //g")
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
