#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/cutandrun
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/cutandrun
    Website: https://nf-co.re/cutandrun
    Slack  : https://nfcore.slack.com/channels/cutandrun
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.fasta) {
    params.bowtie2   = WorkflowMain.getGenomeAttribute(params, 'bowtie2')
} else {
    params.bowtie2   = null
}

params.fasta     = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf       = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gene_bed  = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.blacklist = WorkflowMain.getGenomeAttribute(params, 'blacklist')

/*
========================================================================================
    SPIKEIN GENOME PARAMETER VALUES
========================================================================================
*/

if (!params.spikein_fasta) {
    params.spikein_bowtie2 = WorkflowMain.getGenomeAttributeSpikeIn(params, 'bowtie2')
} else {
    params.spikein_bowtie2 = null
}
params.spikein_fasta   = WorkflowMain.getGenomeAttributeSpikeIn(params, 'fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUTANDRUN } from './workflows/cutandrun'

workflow NFCORE_CUTANDRUN {
    /*
     * WORKFLOW: Run main nf-core/cutandrun analysis pipeline
     */
    CUTANDRUN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * WORKFLOW: Execute a single named workflow for the pipeline
 * See: https://github.com/nf-core/rnaseq/issues/619
 */
workflow {
    NFCORE_CUTANDRUN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
