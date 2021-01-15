////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////



// ////////////////////////////////////////////////////
// /* --          VALIDATE INPUTS                 -- */
// ////////////////////////////////////////////////////

// // TODO nf-core: Add all file path parameters for the pipeline to the list below
// // Check input path parameters to see if they exist
// checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }

// ////////////////////////////////////////////////////
// /* --          CONFIG FILES                    -- */
// ////////////////////////////////////////////////////

// ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// ////////////////////////////////////////////////////
// /* --       IMPORT MODULES / SUBWORKFLOWS      -- */
// ////////////////////////////////////////////////////

// // Don't overwrite global params.modules, create a copy instead and use that within the main script.
// def modules = params.modules.clone()

// def multiqc_options   = modules['multiqc']
// multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// // Local: Modules
// include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions' addParams( options: [publish_files : ['csv':'']] )

// // Local: Sub-workflows
// include { INPUT_CHECK           } from './modules/local/subworkflow/input_check'       addParams( options: [:]                          )

// // nf-core/modules: Modules
// include { FASTQC                } from './modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc']            )
// include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options              )

// ////////////////////////////////////////////////////
// /* --           RUN MAIN WORKFLOW              -- */
// ////////////////////////////////////////////////////

// // Info required for completion email and summary
// def multiqc_report = []

// workflow {

//     ch_software_versions = Channel.empty()

//     /*
//      * SUBWORKFLOW: Read in samplesheet, validate and stage input files
//      */
//     INPUT_CHECK ( 
//         ch_input
//     )

//     /*
//      * MODULE: Run FastQC
//      */
//     FASTQC (
//         INPUT_CHECK.out.reads
//     )
//     ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    

//     /*
//      * MODULE: Pipeline reporting
//      */
//     GET_SOFTWARE_VERSIONS ( 
//         ch_software_versions.map { it }.collect()
//     )

//     /*
//      * MultiQC
//      */  
//     if (!params.skip_multiqc) {
//         workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
//         ch_workflow_summary = Channel.value(workflow_summary)

//         ch_multiqc_files = Channel.empty()
//         ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
//         ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
//         ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//         ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
//         ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        
//         MULTIQC (
//             ch_multiqc_files.collect()
//         )
//         multiqc_report       = MULTIQC.out.report.toList()
//         ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
//     }
// }