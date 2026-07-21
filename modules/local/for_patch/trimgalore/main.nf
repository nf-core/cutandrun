process TRIMGALORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::trim-galore=0.6.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.6--0' :
        'biocontainers/trim-galore:0.6.6--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*trimmed.fastq.gz")  , emit: reads
    tuple val(meta), path("*report.txt")        , emit: log
    path  "versions.yml"                        , emit: versions

    tuple val(meta), path("*.html"), emit: html, optional: true
    tuple val(meta), path("*.zip") , emit: zip, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Clipping presets have to be evaluated in the context of SE/PE
    def c_r1   = params.clip_r1 > 0             ? "--clip_r1 ${params.clip_r1}"                         : ''
    def c_r2   = params.clip_r2 > 0             ? "--clip_r2 ${params.clip_r2}"                         : ''
    def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def suffix = task.ext.suffix ? "${task.ext.suffix}.trimmed" : ".trimmed"
    def prefix_1 = "${meta.id}_1${suffix}"
    def prefix_2 = "${meta.id}_2${suffix}"
    def prefix = "${meta.id}${suffix}"

    def args = task.ext.args ?: ''

    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trim_galore \\
            $args \\
            --cores $cores \\
            --gzip \\
            $c_r1 \\
            $tpc_r1 \\
            ${prefix}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim_galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${meta.id}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.id}_1.fastq.gz
        [ ! -f  ${meta.id}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.id}_2.fastq.gz
        trim_galore \\
            $args \\
            --cores $cores \\
            --paired \\
            --gzip \\
            $c_r1 \\
            $c_r2 \\
            $tpc_r1 \\
            $tpc_r2 \\
            ${meta.id}_1.fastq.gz \\
            ${meta.id}_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trim_galore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS

        mv ${meta.id}_1_val_1.fq.gz ${prefix_1}.fastq.gz
        mv ${meta.id}_2_val_2.fq.gz ${prefix_2}.fastq.gz

        [ ! -f  ${meta.id}_1_val_1_fastqc.html ] || mv ${meta.id}_1_val_1_fastqc.html ${meta.id}_1${suffix}_fastqc.html
        [ ! -f  ${meta.id}_2_val_2_fastqc.html ] || mv ${meta.id}_2_val_2_fastqc.html ${meta.id}_2${suffix}_fastqc.html

        [ ! -f  ${meta.id}_1_val_1_fastqc.zip ] || mv ${meta.id}_1_val_1_fastqc.zip ${meta.id}_1${suffix}_fastqc.zip
        [ ! -f  ${meta.id}_2_val_2_fastqc.zip ] || mv ${meta.id}_2_val_2_fastqc.zip ${meta.id}_2${suffix}_fastqc.zip
        """
    }
}
