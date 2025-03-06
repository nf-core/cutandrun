process PEAK_FRIP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(peaks_bed), path(fragments_bed), path(flagstat)
    path  frip_score_header
    val   min_frip_overlap

    output:
    tuple val(meta), path("*mqc.tsv"), emit: frip_mqc
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    READS_IN_PEAKS=\$(bedtools intersect -a ${fragments_bed} -b ${peaks_bed} -bed -c -f $min_frip_overlap |  awk -F '\t' '{sum += \$NF} END {print sum * 2}')
    grep -m 1 'mapped (' ${flagstat} | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "Peak FRiP Score", a/\$1}' | cat $frip_score_header - > ${prefix}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
