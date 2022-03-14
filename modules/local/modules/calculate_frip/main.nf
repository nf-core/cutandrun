process CALCULATE_FRIP {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::deeptools=3.5.* bioconda::pysam=0.17.*" : null)
    container "luslab/cutandrun-dev-frip:latest"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path('frips.csv'), emit: frips
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    frip.py \\
        --bams "*.bam" \\
        --peaks "*.bed" \\
        --threads ${task.cpus} \\
        --outpath .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """
}
