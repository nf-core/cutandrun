// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEEPTOOLS_BAMPEFRAGMENTSIZE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::deeptools=3.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/deeptools:3.5.0--py_0"
    } else {
        container "quay.io/biocontainers/deeptools:3.5.0--py_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path  blacklist

    output:
    tuple val(meta), path("*_summary.csv"), emit: summary_csv
    tuple val(meta), path("*_raw.csv")    , emit: raw_csv
    tuple val(meta), path("*_log.txt")    , emit: log
    path  "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bamPEFragmentSize \\
        -b $bam \\
        --outRawFragmentLengths ${prefix}_raw.csv \\
        --table ${prefix}_summary.csv \\
        -bl $blacklist \\
        -p $task.cpus \\
        $options.args \\
        > ${prefix}_log.txt
    bamPEFragmentSize --version | sed -e "s/bamPEFragmentSize //g" > ${software}.version.txt
    """
}