process IGV_SESSION {
    tag "igv"
    label 'process_min'

    conda "conda-forge::python=3.8.3"
    container "biocontainers/python:3.8.3"

    input:
    path genome
    path genome_index
    tuple val(meta), path(gtf_bed), path(gtf_bed_index)
    path beds
    path secondary_beds
    path bigwig
    val sort_by_groups

    output:
    path '*.{txt,xml,bed,bigWig,fa,fai,fna,gtf,gff,narrowPeak,broadPeak,gz,tbi,bedGraph}', emit: session, includeInputs: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def colour_pallete = ['38,70,83', '231,111,81', '42,157,143', '244,162,97', '233,196,106']

    def output = ''
    def colours = [:]
    def colour_pos = 0
    def file_list = []

    if (sort_by_groups) {
        file_list = beds.collect { it.toString() }.sort()
        file_list += secondary_beds.collect { it.toString() }.sort()
        file_list += bigwig.collect { it.toString() }.sort()
    }
    else {
        file_list = (bigwig + secondary_beds + beds).collect { it.toString() }.sort()
    }

    file_list.each { file ->
        def group = file.split('_R')[0]
        if (!colours.containsKey(group)) {
            colours[group] = colour_pallete[colour_pos % colour_pallete.size()]
            colour_pos += 1
        }
        output += "${file}\t${colours[group]}\n"
    }

    output = output.trim()
    """
    echo "${output}" > exp_files.txt
    find -L * -iname "*.gtf" -exec echo -e {}"\\t0,48,73" \\; > gtf.igv.txt
    find -L * -iname "*.gff" -exec echo -e {}"\\t0,48,73" \\; > gff.igv.txt
    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ${genome} ${gtf_bed} --path_prefix './'

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
END_VERSIONS
    """
}
