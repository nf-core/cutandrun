#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/cutandrun
========================================================================================
 nf-core/cutandrun Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/cutandrun
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/cutandrun --input samplesheet.csv -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.fasta        = Checks.get_genome_attribute(params, 'fasta')
params.gtf          = Checks.get_genome_attribute(params, 'gtf')
params.blacklist    = Checks.get_genome_attribute(params, 'blacklist')
params.gene_bed     = Checks.get_genome_attribute(params, 'bed12')

log.info 'test-' + params.blacklist

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)


/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow NFCORE_CUTANDRUN {
    /*
     * WORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
     */
    if (params.public_data_ids) {
        include { SRA_DOWNLOAD } from './workflows/sra_download'
        SRA_DOWNLOAD ()
    
    /*
     * WORKFLOW: Run main nf-core/cutandrun analysis pipeline
     */
    } else {
        include { CUTANDRUN } from './workflows/cutandrun' addParams( summary_params: summary_params )
        CUTANDRUN ()
    }
}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

/*
 * WORKFLOW: Execute a single named workflow for the pipeline
 * See: https://github.com/nf-core/rnaseq/issues/619
 */
workflow {
    NFCORE_CUTANDRUN ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
