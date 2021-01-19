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

// Genome
def publish_genome_options = params.save_reference ? [publish_dir: 'genome/target']               : [publish_files: false]
def spikein_genome_options = params.save_reference ? [publish_dir: 'genome/spikein']              : [publish_files: false]
def bowtie2_index_options = params.save_reference ? [publish_dir: 'genome/target/index']          : [publish_files: false]
def bowtie2_spikein_index_options = params.save_reference ? [publish_dir: 'genome/spikein/index'] : [publish_files: false]

// QC
def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

// Trimming
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
if (params.save_trimmed) { trimgalore_options.publish_files.put('fq.gz','') }

// Alignment
def prepareToolIndices  = ['bowtie2']

def bowtie2_align_options          = modules['bowtie2_align']
def bowtie2_spikein_align_options  = modules['bowtie2_spikein_align']
if (params.save_unaligned)         { bowtie2_align_options.publish_files.put('.gz','') }
if (!params.save_spikein_aligned)  { bowtie2_spikein_align_options['publish_files'] = false }

def samtools_sort_options         = modules['samtools_sort']
def samtools_spikein_sort_options = modules['samtools_spikein_sort']
// if (['star_salmon','hisat2'].contains(params.aligner)) {
//     if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// } else {
//     if (params.save_align_intermeds || params.skip_markduplicates) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// }

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
include { PREPARE_GENOME  } from './modules/local/subworkflow/prepare_genome'   addParams( genome_options: publish_genome_options,
                                                                                           spikein_genome_options: spikein_genome_options,
                                                                                           bt2_index_options: bowtie2_index_options,
                                                                                           bt2_spikein_index_options: bowtie2_spikein_index_options,
                                                                                           spikein_fasta: spikein_fasta )
include { ALIGN_BOWTIE2      } from './modules/local/subworkflow/align_bowtie2' addParams( align_options: bowtie2_align_options, 
                                                                                           spikein_align_options: bowtie2_spikein_align_options, 
                                                                                           samtools_options: samtools_sort_options,
                                                                                           samtools_spikein_options: samtools_spikein_sort_options )

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
    ch_trimmed_reads     = FASTQC_TRIMGALORE.out.reads
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    /*
     * SUBWORKFLOW: Alignment to target and spikein genome using botwtie2
     */
    ch_orig_bam                   = Channel.empty()
    ch_orig_spikein_bam           = Channel.empty()
    ch_bowtie2_log                = Channel.empty()
    ch_bowtie2_spikein_log        = Channel.empty()
    ch_samtools_bam               = Channel.empty()
    ch_samtools_bai               = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_samtools_spikein_bam       = Channel.empty()
    ch_samtools_spikein_bai       = Channel.empty()
    ch_samtools_spikein_stats     = Channel.empty()
    ch_samtools_spikein_flagstat  = Channel.empty()
    ch_samtools_spikein_idxstats  = Channel.empty()
    if (params.aligner == 'bowtie2') {
        ALIGN_BOWTIE2 (
            ch_trimmed_reads,
            PREPARE_GENOME.out.bowtie2_index,
            PREPARE_GENOME.out.bowtie2_spikein_index
        )
        ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
        // ch_orig_bam                   = ALIGN_BOWTIE2.out.orig_bam
        // ch_orig_spikein_bam           = ALIGN_BOWTIE2.out.orig_spikein_bam
        // ch_bowtie2_log                = ALIGN_BOWTIE2.out.bowtie2_log
        // ch_bowtie2_spikein_log        = ALIGN_BOWTIE2.out.bowtie2_spikein_log
        // ch_samtools_bam               = ALIGN_BOWTIE2.out.bam
        // ch_samtools_bai               = ALIGN_BOWTIE2.out.bai
        // ch_samtools_stats             = ALIGN_BOWTIE2.out.stats
        // ch_samtools_flagstat          = ALIGN_BOWTIE2.out.flagstat
        // ch_samtools_idxstats          = ALIGN_BOWTIE2.out.idxstats
        // ch_samtools_spikein_bam       = ALIGN_BOWTIE2.out.spikein_bam
        // ch_samtools_spikein_bai       = ALIGN_BOWTIE2.out.spikein_bai 
        // ch_samtools_spikein_stats     = ALIGN_BOWTIE2.out.spikein_stats
        // ch_samtools_spikein_flagstat  = ALIGN_BOWTIE2.out.spikein_flagstat
        // ch_samtools_spikein_idxstats  = ALIGN_BOWTIE2.out.spikein_idxstats
    }

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