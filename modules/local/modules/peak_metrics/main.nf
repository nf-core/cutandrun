process PEAK_METRICS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed) 
    tuple val(meta), path(bam) 
    tuple val(meta), path(flagstat) 
    path  frip_score_header 
    val   min_frip_overlap


    output:
    path '*_mqc.tsv'                  , emit: frip_mqc
    tuple val(meta), path('frips.csv'), emit: frips
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    echo "frip" > frip_header.txt
    READS_IN_PEAKS=\$(bedtools intersect -a ${bam} -b ${bed} -bed -c -f $min_frip_overlap |  awk -F '\t' '{sum += \$NF} END {print sum}')
    grep -m 1 'mapped (' ${flagstat} | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${meta.id}", a/\$1}' | cat $frip_score_header - > ${prefix}_peaks.frip_mqc.tsv
    grep -m 1 'mapped (' ${flagstat} | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print a/\$1}' | cat frip_header.txt - > frips.csv   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
