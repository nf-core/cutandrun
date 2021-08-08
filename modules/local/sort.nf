include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options   = [:]
options          = initOptions(params.options)
options.command  = params.options.command ?: ''
options.ext      = params.options.ext ?: ''


process SORT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }


    conda (params.enable_conda ? null : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
        tuple val(meta), path(input)

    output:
        tuple val(meta), path("*.sort.*"), emit: file

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def ext   = options.ext ? "${options.ext}" : "txt"

    // Option for multiple files to sort
    String input_files = input.join(" ")

    """
    sort $options.args $input_files > ${prefix}.sort.${ext}
    """
}
