// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './common/functions'

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_GENOMECOV_SCALE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    tuple val(meta), path(bam), val(scale)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path  "versions.yml"              , emit: versions
    
    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    bedtools genomecov -ibam $bam $options.args -bg -scale $scale \\
    | bedtools sort > ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
