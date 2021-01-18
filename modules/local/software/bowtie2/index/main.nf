include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE2_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1"
    } else {
        container "quay.io/biocontainers/bowtie2:2.4.2--py37h8270d21_0"
    }

    input:
    path fasta

    output:
    path "*.bt2*"        , emit: index
    path "*.log"         , emit: log
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    bowtie2-build ${options.args} --threads ${task.cpus} $fasta ${fasta.simpleName} > ${fasta.simpleName}.summary.log
    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
