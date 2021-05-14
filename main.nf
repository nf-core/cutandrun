#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/cutandrun
========================================================================================
    Github : https://github.com/nf-core/cutandrun 
    Website: https://nf-co.re/cutandrun
    Slack  : https://nfcore.slack.com/channels/cutandrun
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta     = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf       = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gene_bed  = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.blacklist = WorkflowMain.getGenomeAttribute(params, 'blacklist')
params.bowtie2   = WorkflowMain.getGenomeAttribute(params, 'bowtie2')

/*
========================================================================================
    SPIKEIN GENOME PARAMETER VALUES
========================================================================================
*/

params.spikein_fasta   = WorkflowMain.getGenomeAttributeSpikeIn(params, 'fasta')
params.spikein_bowtie2 = WorkflowMain.getGenomeAttributeSpikeIn(params, 'bowtie2')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)


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

/*
========================================================================================
    THE END
========================================================================================
*/
