////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input,
    params.fasta,
    params.bowtie2_index,
    params.spikein_fasta,
    params.spikein_bowtie2_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }

// Resolve spike-in genome
def spikein_fasta = params.spikein_fasta
if(!params.spikein_bowtie2_index && !params.spikein_fasta) {
    if(params.genomes[params.spikein_geome] != null) {
        spikein_fasta = params.genomes[params.spikein_geome].fasta
    }
    else {
        exit 1, 'Invalid spike-in genome specified!'
    }
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
//ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --               CONFIG FILES               -- */
////////////////////////////////////////////////////



////////////////////////////////////////////////////
/* --     INIALISE PARAMETERS AND OPTIONS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// QC
def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

// Trimming
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

// Alignment
def prepareToolIndices  = ['bowtie2']
def samtools_sort_options = modules['samtools_sort']

// Genome
def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def bowtie2_index_options = modules['bowtie2_index']
if (!params.save_reference) { bowtie2_index_options['publish_files'] = false }
def bowtie2_spikein_index_options = modules['bowtie2_spikein_index']
if (!params.save_reference) { bowtie2_index_options['publish_files'] = false }

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

/*
 * MODULES
 */
include { INPUT_CHECK           } from './modules/local/subworkflow/input_check'       addParams( options: [:] )
include { CAT_FASTQ             } from './modules/local/process/cat_fastq'             addParams( options: cat_fastq_options )
include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions' addParams( options: [publish_files : ['csv':'']] )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
include { PREPARE_GENOME  } from './modules/local/subworkflow/prepare_genome' addParams( genome_options: publish_genome_options, index_options: publish_index_options, bt2_index_options: bowtie2_index_options, bt2_spikein_index_options: bowtie2_spikein_index_options, spikein_fasta: spikein_fasta )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULES
 */


/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { FASTQC_TRIMGALORE } from './modules/nf-core/subworkflow/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow CUTANDRUN {

    // Init
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        prepareToolIndices
    )
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.bowtie2_version.ifEmpty(null))

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_cat_fastq }
    ch_cat_fastq | view

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ ( 
        ch_cat_fastq
    )

    /*
     * SUBWORKFLOW: Read QC, trim adapters and perform post-trim read QC
     */
    FASTQC_TRIMGALORE (
        CAT_FASTQ.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

// workflow.onComplete {
//     Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
//     Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////