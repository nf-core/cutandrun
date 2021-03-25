include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
params.multiqc_label = ''
def options          = initOptions(params.options)

process DESEQ2_QC {
    tag "DESeq2"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.2 conda-forge::r-rcolorbrewer=1.1_2 conda-forge::r-pheatmap=1.0.12 bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel=1.22.0 conda-forge::r-stringr=1.4.0 conda-forge::r-magrittr=2.0.1 bioconda::bioconductor-chromvar=1.10.0 bioconda::bioconductor-genomicranges=1.40.0" : null)
    container "luslab/cutrun-ds2-dev"

    input:
    val groups
    path peak_beds
    path bams

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    //path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path  "*.version.txt"       , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def label_lower = params.multiqc_label.toLowerCase()
    def label_upper = params.multiqc_label.toUpperCase()
    """
    deseq2_diff.r \\
        --groups $groups \\
        --bed $peak_beds \\
        --bam $bams \\
        --cores $task.cpus \\
        $options.args
    """
}