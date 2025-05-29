//
// This file holds several functions specific to the main.nf workflow in the nf-core/cutandrun pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "  https://doi.org/10.5281/zenodo.5653535\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log, args) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)
        // Check that the profile doesn't contain spaces and doesn't end with a trailing comma
        checkProfile(workflow.profile, args, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        println "=== DEBUG getGenomeAttribute ==="
        println "Requested attribute: '${attribute}'"
        println "params.genome: '${params.genome}'"
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            println "params.genomes[ params.genome ]: '${params.genomes[ params.genome ]}'"
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                println "params.genomes[ params.genome ][ attribute ]: '${params.genomes[ params.genome ][ attribute ]}'"
                   return params.genomes[ params.genome ][ attribute ]
            }
            else {
                println "params.genomes[ params.genome ] does not contain attribute: '${attribute}'"
            }
        }
        return null
    }



// public static Object getGenomeAttribute(params, attribute) {
    // log.info "=== DEBUG getGenomeAttribute ==="
    // log.info "Requested attribute: '${attribute}'"
    // log.info "params.genome: '${params.genome}'"
    // 
   // Check if params.genomes exists
    // if (!params.genomes) {
        // log.warn "DEBUG: params.genomes is null or empty"
        // return null
    // }
    // log.info "DEBUG: params.genomes exists and contains ${params.genomes.size()} entries"
    // log.info "DEBUG: Available genomes: ${params.genomes.keySet()}"
    // 
   // Check if params.genome is set
    // if (!params.genome) {
        // log.warn "DEBUG: params.genome is null or empty"
        // return null
    // }
    // log.info "DEBUG: params.genome is set to: '${params.genome}'"
    // 
   // Check if the specified genome exists in params.genomes
    // if (!params.genomes.containsKey(params.genome)) {
        // log.warn "DEBUG: Genome '${params.genome}' not found in available genomes"
        // log.warn "DEBUG: Available genomes are: ${params.genomes.keySet()}"
        // return null
    // }
    //  log.info "DEBUG: Found genome '${params.genome}' in genomes config"
    //  
   // Get the genome config
    // def genomeConfig = params.genomes[params.genome]
    // log.info "DEBUG: Genome config for '${params.genome}': ${genomeConfig}"
    // log.info "DEBUG: Available attributes for '${params.genome}': ${genomeConfig.keySet()}"
    // 
   // Check if the requested attribute exists
    // if (!genomeConfig.containsKey(attribute)) {
        // log.warn "DEBUG: Attribute '${attribute}' not found for genome '${params.genome}'"
        // log.warn "DEBUG: Available attributes are: ${genomeConfig.keySet()}"
        // return null
    // }
    // 
   // Get the attribute value
    // def attributeValue = genomeConfig[attribute]
    // log.info "DEBUG: Found attribute '${attribute}' for genome '${params.genome}'"
    // log.info "DEBUG: Attribute value: '${attributeValue}'"
    // log.info "DEBUG: Attribute value type: ${attributeValue?.getClass()?.getSimpleName()}"
    // log.info "=== END DEBUG getGenomeAttribute ==="
    // 
    // return attributeValue
// }


    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttributeSpikeIn(params, attribute) {
        def val = ''
        if (params.genomes && params.spikein_genome && params.genomes.containsKey(params.spikein_genome)) {
            if (params.genomes[ params.spikein_genome ].containsKey(attribute)) {
                val = params.genomes[ params.spikein_genome ][ attribute ]
            }
        }
        return val
    }

    //
    // Exit pipeline if --profile contains spaces
    //
    private static void checkProfile(profile, args, log) {
        if (profile.endsWith(',')) {
            Nextflow.error "Profile cannot end with a trailing comma. Please remove the comma from the end of the profile string.\nHint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
        if (args[0]) {
            log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${args[0]}` has been detected.\n      Hint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
    }
}
