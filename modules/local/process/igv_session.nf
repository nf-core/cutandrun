include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IGV_SESSION {
    tag "igv"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
      path genome
      path gtf
      path beds
      path bw
    
    output:
      path('*.{txt,xml,bed,bw,fa}', includeInputs:true)

    script:
    """
    find -L * -iname "*.fa" -exec echo -e {}"\\t0,0,178" \\; > fa.igv.txt
    find -L * -iname "*.bed" -exec echo -e {}"\\t0,0,178" \\; > bed.igv.txt
    find -L * -iname "*.bw" -exec echo -e {}"\\t0,0,178" \\; > bw.igv.txt
    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt $genome --path_prefix './'
    """
}