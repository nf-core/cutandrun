include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GENERATE_REPORTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "luslab/cutandrun-dev-reporting:latest"

    input:
    path meta_data
    path meta_data_ctrl
    path raw_fragments
    path bed_fragments
    path seacr_beds
    path frag_len_header_multiqc

    output:
    path '*.pdf'             , emit: pdf
    path '*.csv'             , emit: csv
    path '*.png'             , emit: png
    path '*frag_len_mqc.yaml', emit: frag_len_multiqc
    path  "versions.yml"     , emit: versions

    script:
    def meta_data_resolved = meta_data ? meta_data : meta_data_ctrl

    """
    reporting.py gen_reports \\
        --meta $meta_data_resolved \\
        --meta_ctrl $meta_data_ctrl \\
        --raw_frag "*.frag_len.txt" \\
        --bin_frag "*bin500.awk.bed" \\
        --seacr_bed "*bed.*.bed" \\
        --output . \\
        --log log.txt

    if [ -f "03_03_frag_len_mqc.txt" ]; then
        cat $frag_len_header_multiqc 03_03_frag_len_mqc.txt > frag_len_mqc.yaml
    fi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
    END_VERSIONS
    """

}
