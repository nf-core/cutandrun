//
// This file holds several functions specific to the workflow/cutandrun.nf in the nf-core/cutandrun pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowCutandrun {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)

        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }

        if (params.normalisation_mode == "Spikein" && !params.spikein_fasta) {
            def error_string = "Spike-in fasta file not specified with e.g. '--spikein_fasta genome.fa' or via a detectable config file."
            Nextflow.error(error_string)
        }

        if (!params.gtf) {
            def error_string = "No GTF annotation specified!"
            Nextflow.error(error_string)
        }

        if (params.gtf) {
            if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
                ncbiGenomeWarn(log)
            }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }

        if (params.dt_calc_all_matrix) {
            matrixWarn(log)
        }

        if (!params.genome && !params.mito_name && params.remove_mitochondrial_reads) {
            rmMitoWarn(log)
        }

        if (!params.fasta && !params.mito_name && params.remove_mitochondrial_reads) {
            rmMitoWarn(log)
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
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO nf-core: Optionally add in-text citation tools to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Print a warning if using GRCh38 assembly from igenomes.config
    //
    private static void ncbiGenomeWarn(log) {
        log.warn "=============================================================================\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
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

    private static void varryingReplicateNumbersError(log) {
        def error_string = "===================================================================================\n" +
            "  There are varrying numbers of replicates across experiemental and IgG samples.\n" +
            "  Options:\n" +
            "    - provide a consistent number of replicates across all experiments and control\n" +
            "    - provide any number of experimental replicates along with just one control rep\n" +
            "    - provide any number of experimental replicates and give all control replicates\n" +
            "      the same replicate number, so that they will be merged for downstream analysis\n" +
            "==================================================================================="
        Nextflow.error(error_string)
    }

    private static void varryingReplicateNumbersWarn(log) {
        log.warn "===================================================================================\n" +
            "  The number of replicates for IgG control does not match the number of replicates \n" +
            "  for experimental data. Only the first IgG replicate will be used for SEACR \n" +
            "  peak-caller normalisation and downstream analysis.\n" +
            "==================================================================================="
    }

    //
    // Print a warning if not blacklist detected
    //
    private static void blacklistWarn(log) {
        log.warn "=============================================================================\n" +
            "  No genome blacklist file specified, switching to dummy empty file...\n" +
            "==================================================================================="
    }

    //
    // Print a warning if gen all plots are on
    //
    private static void matrixWarn(log) {
        log.warn "=========================================================================================================\n" +
            "  dt_calc_all_matrix is switched on which will calculate a deeptools matrix for all samples. \n" +
            "  If you have a large sample count, this may affect pipeline performance and result in errors. \n" +
            "  Set this option to false to disable this feature and only calculate deeptools heatmaps for single samples\n" +
            "==============================================================================================================="
    }

    //
    // Print a warning if neither genome or mito_name are specified and user wants to filter mitochodnrial reads
    //
    private static void rmMitoWarn(log) {
        log.warn "=============================================================================\n" +
            "  No genome or mito_name specified but remove_mitochondrial_reads set to true \n" +
            "  remove_mitochondrial_reads is now set to false as the name of mitochondrial reads \n" +
            "  in the .bam files is unknown. \n" +
            "==================================================================================="
    }
}
