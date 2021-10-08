include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALCULATE_FRIP {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3 bioconda::deeptools=3.5.* bioconda::pysam=0.17.*" : null)
    container "luslab/cutandrun-dev-frip:latest"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path('frips.csv'), emit: frips
    path '*.version.txt'              , emit: version

    script:
    """
    frip.py \\
        --bams "*.bam" \\
        --peaks "*.bed" \\
        --threads ${task.cpus} \\
        --outpath .

    python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\" > python.version.txt
    """
}
