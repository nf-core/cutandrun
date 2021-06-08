//
// This file holds several functions specific to the workflow/cutandrun.nf in the nf-core/cutandrun pipeline
//

class WorkflowCutandrun {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)

        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.spikein_fasta) {
            log.error "Spike-in fasta file not specified with e.g. '--spikein_fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.gtf) {
            log.error "No GTF annotation specified!"
            System.exit(1)
        }

        if (params.gtf) {
            if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
                ncbiGenomeWarn(log)
            }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }

    //
    // Print a warning if using GRCh38 assembly from igenomes.config
    //
    private static void ncbiGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    //
    // Print a warning if using a UCSC assembly from igenomes.config
    //
    private static void ucscGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "==================================================================================="
    }

    //
    // IgG workflow channel mapping function
    //

    private static String[][] checkReplicateNumbers(row) {
        def exp_reps = row.last()
        def igg_reps = row[row.size() - 2]
        def current_reps = row[0]
        def unique_exp_reps = exp_reps.unique(false)
        def unique_igg_reps = igg_reps.unique(false)
        def exp_rep_freq = [0] * unique_exp_reps.size()
        def output = row[1]
        def final_output = []
        def all_same = false
        def i_freq = 0

        // check if exp rep numbers are occuring an equal number of times
        for (i=0; i<unique_exp_reps.size(); i++) {
            i_freq = exp_reps.count(unique_exp_reps[i])
            exp_rep_freq[i] = i_freq
        }
        all_same = exp_rep_freq.every{ it ==  exp_rep_freq[0] }

        // check cases and assign if criteria is met
        if ( all_same && (unique_exp_reps.sort() ==  unique_igg_reps.sort()) && (current_reps[0] == current_reps[1]) ) {
            final_output = output
        } else if ( unique_igg_reps.size() == 1 ) {
            final_output = output
        } else if ( all_same && (unique_igg_reps.size() != 1) && (current_reps[1] == unique_igg_reps.min()) ) {
            WorkflowCutandrun.varryingReplicateNumbersWarn(log)
            final_output = output
        } else if ( !all_same && (unique_igg_reps.size() != 1) ) {
            WorkflowCutandrun.varryingReplicateNumbersError(log)
        }

        return final_output
    }

    private static void varryingReplicateNumbersError(log) {
        log.error "=============================================================================\n" +
            "  There are varrying numbers of replicates across experiemental and IgG samples.\n" +
            "  Options:\n" +
            "    - provide a consistent number of replicates across all experiments and control\n" +
            "    - provide any number of experimental replicates along with just one control rep\n" +
            "    - provide any number of experimental replicates and give all control replicates\n" +
            "      the same replicate number, so that they will be merged for downstream analysis\n" +
            "==================================================================================="
        System.exit(1)
    }

    private static void varryingReplicateNumbersWarn(log) {
        log.warn "=============================================================================\n" +
            "  The number of replicates for IgG control does not match the number of replicates \n" +
            "  for experimental data. Only the first IgG replicate will be used for SEACR \n" +
            "  peak-caller normalisation and downstream analysis.\n" +
            "==================================================================================="
    }

}
