include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
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

    output:

    script:
    
}