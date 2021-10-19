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
    path raw_fragments
    path bed_fragments
    path seacr_beds
    path frag_len_header_multiqc

    output:
    path '*.pdf',             emit: pdf
    path '*.csv',             emit: csv
    path '*.png',             emit: png
    path '*frag_len_mqc.yaml', emit: frag_len_multiqc
    path '*.version.txt',     emit: version

    script:  // This script is bundled with the pipeline, in nf-core/cutandrun/bin/
    """
    reporting.py gen_reports \\
        --meta $meta_data \\
        --raw_frag "*.frag_len.txt" \\
        --bin_frag "*bin500.awk.bed" \\
        --seacr_bed "*bed.*.bed" \\
        --output . \\
        --log log.txt

    if [ -f "03_03_frag_len_mqc.txt" ]; then
        cat $frag_len_header_multiqc 03_03_frag_len_mqc.txt > frag_len_mqc.yaml
    fi

    python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\" > python.version.txt
    """

}
