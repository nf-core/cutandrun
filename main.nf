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

log.info Headers.nf_core(workflow, params.monochrome_logs)

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
        include { CUTANDRUN } from './workflows/cutandrun'
        CUTANDRUN ()
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/cutandrun]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/cutandrun]${c_red} Pipeline completed with errors${c_reset}-"
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
