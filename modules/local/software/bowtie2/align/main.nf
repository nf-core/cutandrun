include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.save_unaligned = false
def options    = initOptions(params.options)

process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // bowtie2=2.4.2, samtools=1.11, pigz=2.3.4
    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2 bioconda::bowtie2=2.4.2 conda-forge::pigz=2.3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam"),           emit: bam
    tuple val(meta), path("*.unmapped.1.gz"), path("*.unmapped.2.gz"), optional: true, emit: unmapped_sam
    tuple val(meta), path("*summary.txt"),    emit: log
    path  "*.version.txt",                    emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args = "--threads $task.cpus --no-unal"
    def inputs     = ''
    def unmapped  = ''

    read_list = reads.collect{it.toString()}
    if(read_list.size > 1){
        inputs = '-1 ' + reads[0] + ' -2 ' + reads[1]
        if(params.save_unaligned) { unmapped = " --un-conc-gz ${prefix}.unmapped.gz" }
    }
    else {
        inputs = '-U ' + reads[0]
        if(params.save_unaligned) { unmapped = " --un-gz ${prefix}.unmapped.gz" }
    }

    args += unmapped
    if(options.args != '') { args += ' ' + options.args }

    """
    bowtie2 -x ${index[0].simpleName} $args $inputs 2>${prefix}_summary.txt \\
    | samtools view $options.args2 -@ $task.cpus -bS -o ${prefix}.bam -

    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
