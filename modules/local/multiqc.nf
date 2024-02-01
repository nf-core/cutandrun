process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

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
    path ('peak_metrics/peak_count_consensus/*')
    path ('peak_metrics/peak_reprod_perc/*')
    path ('frag_len/*')
    path ('linear_duplicates/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $args $custom_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
