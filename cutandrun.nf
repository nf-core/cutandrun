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
    params.gtf,
    params.blacklist,
    params.bowtie2_index,
    params.spikein_fasta,
    params.spikein_bowtie2_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
if (params.gtf) { ch_gtf = file(params.gtf)   } else { exit 1, 'Genome GTF file not specified!' }
if (params.blacklist) { ch_blacklist = file(params.blacklist) } else { exit 1, 'Genome blacklist file not specified!' }

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

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --                  ASSETS                  -- */
////////////////////////////////////////////////////

ch_dummy_file = file("dummy/file", checkIfExists: false)
ch_bt2_to_csv_awk = file("$projectDir/assets/awk/bt2_report_to_csv.awk", checkIfExists: true)
ch_dt_frag_to_csv_awk = file("$projectDir/assets/awk/dt_frag_report_to_csv.awk", checkIfExists: true)

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

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Trimming
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
if (params.save_trimmed) { trimgalore_options.publish_files.put('fq.gz','') }

// Alignment dedup and filtering
def prepareToolIndices  = ['bowtie2']

def bowtie2_align_options          = modules['bowtie2_align']
def bowtie2_spikein_align_options  = modules['bowtie2_spikein_align']
if (params.save_unaligned)         { bowtie2_align_options.publish_files.put('.gz','') }
if (params.save_unaligned)         { bowtie2_spikein_align_options.publish_files.put('.gz','') }
if (params.publish_align_intermed) { bowtie2_align_options.publish_files.put('.bam','') }
if (params.publish_align_intermed) { bowtie2_spikein_align_options.publish_files.put('.bam','') }

def samtools_sort_options         = modules['samtools_sort']
def samtools_spikein_sort_options = modules['samtools_spikein_sort']
if (params.publish_align_intermed || params.skip_markduplicates) {
    samtools_sort_options.publish_files.put('bam','')
    samtools_sort_options.publish_files.put('bai','')
    samtools_spikein_sort_options.publish_files.put('bam','')
    samtools_spikein_sort_options.publish_files.put('bai','')
}


def samtools_view_options         = modules['samtools_view']
def samtools_qfilter_options         = modules['samtools_qfilter']
if (params.minimum_alignment_q_score > 0) {
    samtools_view_options.args = "-q " + params.minimum_alignment_q_score
}
if (params.publish_align_intermed) {
    samtools_view_options.publish_files = ['bam':'']
    samtools_qfilter_options.publish_files.put('bai','')
}

def picard_markduplicates_options         = modules['picard_markduplicates']
def picard_markduplicates_samtools_options         = modules['picard_markduplicates_samtools']
if (params.publish_align_intermed) {
    picard_markduplicates_options.publish_files.put('bam','')
    picard_markduplicates_samtools_options.publish_files.put('bai','')
}
dedup_control_only = true
if(params.dedup_target_reads) { dedup_control_only = false }

// Meta annotation
def awk_bt2_options = modules['awk_bt2']
def awk_bt2_spikein_options = modules['awk_bt2_spikein']
def awk_dedup_options = modules['awk_dedup']
def awk_dt_frag_options = modules['awk_dt_frag']

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

/*
 * MODULES
 */
include { INPUT_CHECK              } from './modules/local/subworkflow/input_check'          addParams( options: [:] )
include { CAT_FASTQ                } from './modules/local/process/cat_fastq'                addParams( options: cat_fastq_options )
include { BEDTOOLS_GENOMECOV_SCALE } from './modules/local/process/bedtools_genomecov_scale' addParams( options: modules['bedtools_genomecov_bedgraph'] )
include { SEACR_CALLPEAK           } from './modules/local/software/seacr/callpeak/main'     addParams( options: modules['seacr'] )
include { UCSC_BEDCLIP             } from './modules/local/process/ucsc_bedclip'             addParams( options: modules['ucsc_bedclip']  )
include { IGV_SESSION              } from './modules/local/process/igv_session'              addParams( options: modules['igv']  )
include { GET_SOFTWARE_VERSIONS    } from './modules/local/process/get_software_versions'    addParams( options: [publish_files : ['csv':'']] )
include { MULTIQC                            } from './modules/local/process/multiqc'                     addParams( options: multiqc_options )
include { EXPORT_META                            } from './modules/local/process/export_meta'                     addParams( options: modules['export_meta'] )
include { GENERATE_REPORTS                            } from './modules/local/process/generate_reports'                     addParams( options: modules['generate_reports'] )
include { DEEPTOOLS_BAMPEFRAGMENTSIZE } from './modules/local/software/deeptools/bamPEFragmentSize/main' addParams( options: modules['deeptools_fragmentsize'] )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
include { PREPARE_GENOME } from './modules/local/subworkflow/prepare_genome' addParams( genome_options: publish_genome_options,
                                                                                           spikein_genome_options: spikein_genome_options,
                                                                                           bt2_index_options: bowtie2_index_options,
                                                                                           bt2_spikein_index_options: bowtie2_spikein_index_options,
                                                                                           spikein_fasta: spikein_fasta )
include { ALIGN_BOWTIE2 } from './modules/local/subworkflow/align_bowtie2'   addParams( align_options: bowtie2_align_options, 
                                                                                           spikein_align_options: bowtie2_spikein_align_options, 
                                                                                           samtools_options: samtools_sort_options,
                                                                                           samtools_spikein_options: samtools_spikein_sort_options )
include { SAMTOOLS_VIEW_SORT_STATS } from './modules/local/subworkflow/samtools_view_sort_stats' addParams( samtools_options: samtools_qfilter_options, samtools_view_options: samtools_view_options)                                                                                    
include { ANNOTATE_META_AWK as ANNOTATE_BT2_META } from './modules/local/subworkflow/annotate_meta_awk' addParams( options: awk_bt2_options, meta_suffix: '_target', script_mode: true)
include { ANNOTATE_META_AWK as ANNOTATE_BT2_SPIKEIN_META } from './modules/local/subworkflow/annotate_meta_awk' addParams( options: awk_bt2_spikein_options, meta_suffix: '_spikein', script_mode: true)
include { ANNOTATE_META_AWK as ANNOTATE_DEDUP_META } from './modules/local/subworkflow/annotate_meta_awk' addParams( options: awk_dedup_options, meta_suffix: '', meta_prefix: 'dedup_', script_mode: false)
include { ANNOTATE_META_AWK as ANNOTATE_DT_FRAG_META } from './modules/local/subworkflow/annotate_meta_awk' addParams( options: awk_dt_frag_options, meta_suffix: '', meta_prefix: '', script_mode: true)                                                                    

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULES
 */
include { UCSC_BEDGRAPHTOBIGWIG } from './modules/nf-core/software/ucsc/bedgraphtobigwig/main' addParams( options: modules['ucsc_bedgraphtobigwig'] )

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { FASTQC_TRIMGALORE                                        } from './modules/nf-core/subworkflow/fastqc_trimgalore'      addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )
include { MARK_DUPLICATES_PICARD                                   } from './modules/nf-core/subworkflow/mark_duplicates_picard' addParams( markduplicates_options: picard_markduplicates_options, samtools_options: picard_markduplicates_samtools_options, control_only: false )
include { MARK_DUPLICATES_PICARD as DEDUP_PICARD                   } from './modules/nf-core/subworkflow/mark_duplicates_picard' addParams( markduplicates_options: modules['picard_dedup'], samtools_options: modules['picard_dedup_samtools'], control_only: dedup_control_only )

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
        ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
        ch_orig_bam                   = ALIGN_BOWTIE2.out.orig_bam
        ch_orig_spikein_bam           = ALIGN_BOWTIE2.out.orig_spikein_bam
        ch_bowtie2_log                = ALIGN_BOWTIE2.out.bowtie2_log
        ch_bowtie2_spikein_log        = ALIGN_BOWTIE2.out.bowtie2_spikein_log

        ch_samtools_bam               = ALIGN_BOWTIE2.out.bam
        ch_samtools_bai               = ALIGN_BOWTIE2.out.bai
        ch_samtools_stats             = ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat          = ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats          = ALIGN_BOWTIE2.out.idxstats

        ch_samtools_spikein_bam       = ALIGN_BOWTIE2.out.spikein_bam
        ch_samtools_spikein_bai       = ALIGN_BOWTIE2.out.spikein_bai 
        ch_samtools_spikein_stats     = ALIGN_BOWTIE2.out.spikein_stats
        ch_samtools_spikein_flagstat  = ALIGN_BOWTIE2.out.spikein_flagstat
        ch_samtools_spikein_idxstats  = ALIGN_BOWTIE2.out.spikein_idxstats
    }

    /*
     *  SUBWORKFLOW: Filter reads based on quality metrics
     *  http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
     */
    if (params.minimum_alignment_q_score > 0) {
        SAMTOOLS_VIEW_SORT_STATS (
            ch_samtools_bam
        )
        ch_samtools_bam           = SAMTOOLS_VIEW_SORT_STATS.out.bam
        ch_samtools_bai           = SAMTOOLS_VIEW_SORT_STATS.out.bai
        ch_samtools_stats         = SAMTOOLS_VIEW_SORT_STATS.out.stats
        ch_samtools_flagstat      = SAMTOOLS_VIEW_SORT_STATS.out.flagstat
        ch_samtools_idxstats      = SAMTOOLS_VIEW_SORT_STATS.out.idxstats
    }

    /*
     * SUBWORKFLOW: Mark duplicates on all samples
     */
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_samtools_bam
        )
        ch_samtools_bam           = MARK_DUPLICATES_PICARD.out.bam
        ch_samtools_bai           = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Remove duplicates
     */
    ch_dedup_multiqc = Channel.empty()
    if (!params.skip_markduplicates && !params.skip_removeduplicates) {
        DEDUP_PICARD (
            ch_samtools_bam
        )
        ch_samtools_bam           = DEDUP_PICARD.out.bam
        ch_samtools_bai           = DEDUP_PICARD.out.bai
        ch_samtools_stats         = DEDUP_PICARD.out.stats
        ch_samtools_flagstat      = DEDUP_PICARD.out.flagstat
        ch_samtools_idxstats      = DEDUP_PICARD.out.idxstats
        ch_dedup_multiqc          = DEDUP_PICARD.out.metrics
    }

    /*
     * SUBWORKFLOW: Annotate meta data with aligner stats and 
     */
    ANNOTATE_BT2_META( ch_samtools_bam, ch_bowtie2_log, ch_bt2_to_csv_awk)
    ANNOTATE_BT2_SPIKEIN_META( ch_samtools_bam, ch_bowtie2_spikein_log, ch_bt2_to_csv_awk)

    /*
     * CHANNEL: Combine merge spikein meta data with main data stream
     */
    ANNOTATE_BT2_SPIKEIN_META.out.output
        .map { row -> [row[0].id, row[0] ].flatten()}
        .set { ch_spikein_bt2_meta }

    ANNOTATE_BT2_META.out.output
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_spikein_bt2_meta )
        .map { row -> [ row[1] << row[3], row[2] ] }
        .set { ch_combined_meta }

    /*
     * CHANNEL: Calculate scale factor for each sample and join to main data flow
     */
    ANNOTATE_BT2_SPIKEIN_META.out.output
        .map { row -> [ row[0].id, params.normalisation_c / (row[0].find{ it.key == "bt2_total_aligned_spikein" }?.value.toInteger()) ] }
        .set { ch_scale_factor }

    ch_combined_meta
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_scale_factor )
        .map { row -> row[1..(row.size() - 1)] }
        .map { row -> 
            row[0].put('scale_factor', row[2])
            [ row[0], row[1], row[2] ] }
        .set { ch_samtools_bam_scale }

    //Create channel without scale as seperate value
    ch_samtools_bam_scale
        .map { row -> [ row[0], row[1] ] }
        .set { ch_samtools_bam_meta }

    if(!params.skip_coverage) {
        /*
        * MODULE: Convert to bedgraph
        */
        BEDTOOLS_GENOMECOV_SCALE (
            ch_samtools_bam_scale
        )
        ch_bedtools_bedgraph = BEDTOOLS_GENOMECOV_SCALE.out.bedgraph
        ch_software_versions = ch_software_versions.mix(BEDTOOLS_GENOMECOV_SCALE.out.version.first().ifEmpty(null))

        /*
         * CHANNEL: Separate bedgraphs into target/control pairings for each replicate
         */
         ch_bedtools_bedgraph.branch { it ->
            target: it[0].group != 'igg'
            control: it[0].group == 'igg'
        }
        .set { ch_bedgraph_split }

        ch_bedgraph_split.target
            .combine(ch_bedgraph_split.control)
            .filter { row -> row[0].replicate == row[2].replicate }
            .map { row -> [ row[0], row[1], row[3] ] }
            .set { ch_bedgraph_combined }

        /*
        * MODULE: Call peaks
        */
        SEACR_CALLPEAK (
            ch_bedgraph_combined
        )
        ch_seacr_bed = SEACR_CALLPEAK.out.bed
        ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK.out.version.first().ifEmpty(null))

        /*
        * MODULE: Clip off-chromosome peaks
        */
        UCSC_BEDCLIP (
            ch_bedtools_bedgraph,
            PREPARE_GENOME.out.chrom_sizes
        )

        /*
        * MODULE: Convert to bigwig
        */
        UCSC_BEDGRAPHTOBIGWIG (
            UCSC_BEDCLIP.out.bedgraph,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_software_versions = ch_software_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.version.first().ifEmpty(null))

        /*
        * MODULE: Create igv session
        */
        IGV_SESSION (
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf,
            SEACR_CALLPEAK.out.bed.collect{it[1]}.ifEmpty([]),
            UCSC_BEDGRAPHTOBIGWIG.out.bigwig.collect{it[1]}.ifEmpty([])
        )
    }

    /*
     * MODULE: Collect software versions used in pipeline
     */
    // GET_SOFTWARE_VERSIONS ( 
    //     ch_software_versions.map { it }.collect()
    // )

    /*
     * MODULE: Multiqc
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

    /*
     * MODULE: Reporting
     */
    if (!params.skip_reporting) {
        ANNOTATE_DEDUP_META(ch_samtools_bam_meta, ch_markduplicates_multiqc, ch_dummy_file.collect())
        //ANNOTATE_DEDUP_META.out.output | view

        ch_samtools_bam
            .map { row -> [row[0].id, row[0], row[1] ] }
            .set { ch_samtools_bam_id }

        ch_samtools_bai
            .map { row -> [row[0].id, row[0], row[1] ] }
            .set { ch_samtools_bai_id }

        ch_samtools_bam_id
            .join( ch_samtools_bai_id )
            .map { row -> [row[1], row[2], row[4] ] }
            .set { ch_samtools_bam_bai }
        //ch_samtools_bam_bai | view
        
        DEEPTOOLS_BAMPEFRAGMENTSIZE(ch_samtools_bam_bai, ch_blacklist)
        //DEEPTOOLS_BAMPEFRAGMENTSIZE.out.summary_csv | view

        ANNOTATE_DT_FRAG_META( ANNOTATE_DEDUP_META.out.output, DEEPTOOLS_BAMPEFRAGMENTSIZE.out.summary_csv, ch_dt_frag_to_csv_awk)
        //ANNOTATE_DT_FRAG_META.out.output | view

        EXPORT_META (
            ANNOTATE_DEDUP_META.out.output.collect{it[0]}.ifEmpty(['{{NO-DATA}}'])
            //ch_samtools_bam_scale.collect{it[0]}.ifEmpty(['{{NO-DATA}}'])
        )

        GENERATE_REPORTS(EXPORT_META.out.csv, DEEPTOOLS_BAMPEFRAGMENTSIZE.out.raw_csv.collect{it[1]})
    }
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