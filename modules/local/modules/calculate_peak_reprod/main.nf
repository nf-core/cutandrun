include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../common/functions'

params.options = [:]
options        = initOptions(params.options)

process CALCULATE_PEAK_REPROD {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::dask=2021.9.1 conda-forge::pandas=1.3.3" : null)
    container "luslab/cutandrun-dev-peakrepo:latest"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path('peak_repro.csv'), emit: csv
    path  "versions.yml"                   , emit: versions

    script:
    """
    peak_reproducability.py \\
        --intersect $bed \\
        --threads ${task.cpus} \\
        --outpath .

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """
}
