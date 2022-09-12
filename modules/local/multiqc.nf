process MULTIQC {
    label 'process_ultralow'

    conda (params.enable_conda ? "bioconda::multiqc=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path ('software_versions/*')
    path ('software_versions/*')
    path ('workflow_summary/*')
    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')
    path ('bowtie2/*')
    path ('bowtie2_spikein/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('picard/markduplicates/*')
    path ('preseq/*')
    path ('deeptools/*')
    path ('deeptools/*')
    path ('deeptools/*')
    path ('peak_metrics/peak_count/*')
    path ('peak_metrics/peak_frip/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $args $custom_config .
    """
}
