include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

colour_pallete = ['38,70,83', '231,111,81', '42,157,143', '244,162,97', '233,196,106']

process IGV_SESSION {
    tag "igv"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
      path genome
      path gtf
      path beds
      path bigwig
    
    output:
      path('*.{txt,xml,bed,bigWig,fa,gtf}', includeInputs:true)

    script:
    output = ''
    colours = [:]
    colour_pos = 0

    file_list = beds.collect{it.toString()}.sort()
    file_list += bigwig.collect{it.toString()}.sort()
    for(file in file_list){
        file_split = file.split('_R')
        group = file_split[0]
        if(!colours.containsKey(group)) {
            colours.put(group, colour_pallete[colour_pos])

            colour_pos++
            if(colour_pos == colour_pallete.size) {
                colour_pos = 0
            }
        }

        line = file + "\t" + colours[group] + "\n"
        output += line
    }
    output = output.trim()

    """
    echo "$output" > exp_files.txt
    find -L * -iname "*.gtf" -exec echo -e {}"\\t0,48,73" \\; > gtf.igv.txt
    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt $genome --path_prefix './'
    """
}

// find -L * -iname "*.gtf" -exec echo -e {}"\\t0,0,178" \\; > gtf.igv.txt
// find -L * -iname "*.bed" -exec echo -e {}"\\t0,0,178" \\; > bed.igv.txt
// find -L * -iname "*.bigWig" -exec echo -e {}"\\t0,0,178" \\; > bigwig.igv.txt
// cat *.txt > igv_files.txt