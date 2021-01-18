include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
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
    tuple val(meta), path("*.bam"), emit: bam
    //path  "*.version.txt"         , emit: version
    //tuple val(meta), path("${prefix}${opts.unmapped_suffix}.1.fastq.gz"), path("${prefix}${opts.unmapped_suffix}.2.fastq.gz"), optional: true, emit: unmapped_fq_pe
    //tuple val(meta), path("${prefix}${opts.unmapped_suffix}.fastq.gz"), optional: true, emit: unmapped_fq_s
    //tuple val(meta), path("${summary_name}.txt"), emit: report_meta

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputs     = ''
    def unpaired  = ''

    read_list = reads.collect{it.toString()}
    if(read_list.size > 1){
        inputs = '-1 ' + reads[0] + ' -2 ' + reads[1]
        unpaired = '--un-gz'
    }
    else {
        inputs = '-U ' + reads[0]
        unpaired = '--un-conc-gz'
    }

    """
    bowtie2 \\
    $options.args \\
    -p $task.cpus \\
    --no-unal \\
    $unpaired \\
    -x ${index[0].simpleName} \\
    $inputs \\
    | samtools view $options.args2 -@ $task.cpus -bS -o ${prefix}.bam -
    """