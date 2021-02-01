include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process EXPORT_META {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val meta
    
    output:
    path "meta_table.csv", emit: csv
    
    script:
    arr_str = meta[0].keySet().join(",")

    for (int i = 0; i < meta.size(); i++) {
        sample_str = meta[i].values().join(",")
        arr_str =  arr_str + "\n" + sample_str
    }

    """
    echo "${arr_str}" > meta_table.csv
    """
}