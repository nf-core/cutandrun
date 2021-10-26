include { initOptions; saveFiles; getSoftwareName } from './common/functions'

params.options = [:]
options        = initOptions(params.options)

process EXPORT_META {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? null : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val meta
    val table_name

    output:
    path "*.csv", emit: csv

    script:
    def header = [:]

    // Find the header key set
    for (int i = 0; i < meta.size(); i++) {
        meta[i].each {
            entry ->
                if(!header.containsKey(entry.key)) {
                    header.put(entry.key, null)
                }
        }
    }

    // Init output string
    arr_str = header.keySet().join(",")

    // Map the values and write row
    for (int i = 0; i < meta.size(); i++) {
        header.each {
            entry ->
                entry.value = null
        }

        meta[i].each {
            entry ->
                header[entry.key] = entry.value
        }
        sample_str = header.values().join(",")
        arr_str =  arr_str + "\n" + sample_str
    }

    """
    echo "$arr_str" > ${table_name}.csv
    """
}
