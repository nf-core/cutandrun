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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta           = getGenomeAttribute('fasta')
params.bowtie2         = getGenomeAttribute('bowtie2')
params.gtf             = getGenomeAttribute('gtf')
params.gene_bed        = getGenomeAttribute('bed12')
params.blacklist       = getGenomeAttribute('blacklist')
params.spikein_fasta   = getGenomeAttributeSpikeIn('fasta')
params.spikein_bowtie2 = getGenomeAttributeSpikeIn('bowtie2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUTANDRUN               } from './workflows/cutandrun'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cutandrun_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cutandrun_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_CUTANDRUN {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    def caller_list = ['seacr', 'macs2']
    callers = params.peakcaller ? params.peakcaller.split(',').collect { it.trim().toLowerCase() } : ['seacr']
    if ((caller_list + callers).unique().size() != caller_list.size()) {
        error("Invalid variant caller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}")
    }

    //
    // WORKFLOW: Run pipeline
    //
    CUTANDRUN(
        samplesheet,
        params.blacklist ? (Channel.from(file(params.blacklist, checkIfExists: true))) : Channel.empty(),
        file("${projectDir}/bin/bt2_report_to_csv.awk", checkIfExists: true),
        file("${projectDir}/assets/dummy_file.txt", checkIfExists: true),
        ["bowtie2"],
        callers,
    )

    emit:
    multiqc_report = CUTANDRUN.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CUTANDRUN(
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CUTANDRUN.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

def getGenomeAttributeSpikeIn(attribute) {
    if (params.genomes && params.spikein_genome && params.genomes.containsKey(params.spikein_genome)) {
        if (params.genomes[params.spikein_genome].containsKey(attribute)) {
            return params.genomes[params.spikein_genome][attribute]
        }
    }
    return null
}
