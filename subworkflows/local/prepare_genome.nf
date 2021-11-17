/*
 * Uncompress and prepare reference genome files
*/

params.genome_options            = [:]
params.spikein_genome_options    = [:]
params.bt2_index_options         = [:]
params.bt2_spikein_index_options = [:]

include { GUNZIP as GUNZIP_FASTA                     } from "../../modules/nf-core/modules/gunzip/main.nf"     addParams( options: params.genome_options            )
include { GUNZIP as GUNZIP_SPIKEIN_FASTA             } from "../../modules/nf-core/modules/gunzip/main.nf"     addParams( options: params.spikein_genome_options    )
include { GUNZIP as GUNZIP_GTF                       } from "../../modules/nf-core/modules/gunzip/main.nf"     addParams( options: params.genome_options            )
include { GUNZIP as GUNZIP_BED                       } from "../../modules/nf-core/modules/gunzip/main.nf"     addParams( options: params.genome_options            )
include { CUSTOM_GETCHROMSIZES                            } from "../../modules/nf-core/modules/custom/getchromsizes/main.nf"              addParams( options: params.genome_options            )
include { CUSTOM_GETCHROMSIZES as GET_SPIKEIN_CHROM_SIZES } from "../../modules/nf-core/modules/custom/getchromsizes/main.nf"               addParams( options: params.spikein_genome_options    )
include { UNTAR as UNTAR_BT2_INDEX                   } from "../../modules/nf-core/modules/untar/main.nf"      addParams( options: params.bt2_index_options         )
include { UNTAR as UNTAR_SPIKEIN_BT2_INDEX           } from "../../modules/nf-core/modules/untar/main.nf"      addParams( options: params.bt2_spikein_index_options )
include { BOWTIE2_BUILD                              } from "../../modules/nf-core/modules/bowtie2/build/main" addParams( options: params.bt2_index_options         )
include { BOWTIE2_BUILD as BOWTIE2_SPIKEIN_BUILD     } from "../../modules/nf-core/modules/bowtie2/build/main" addParams( options: params.bt2_spikein_index_options )

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:
    ch_versions = Channel.empty()

    /*
    * Uncompress genome fasta file if required
    */
    if (params.fasta.endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( params.fasta ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
    * Uncompress spike-in genome fasta file if required
    */
    if (params.spikein_fasta.endsWith(".gz")) {
        ch_spikein_fasta = GUNZIP_SPIKEIN_FASTA ( params.spikein_fasta ).gunzip
        ch_versions      = ch_versions.mix(GUNZIP_SPIKEIN_FASTA.out.versions)
    } else {
        ch_spikein_fasta = file(params.spikein_fasta)
    }

    /*
    * Uncompress GTF annotation file
    */
    ch_gtf = Channel.empty()
    if (params.gtf.endsWith(".gz")) {
        ch_gtf      = GUNZIP_GTF ( params.gtf ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = file(params.gtf)
    }

    /*
    * Uncompress BED annotation file
    */
    ch_gene_bed = Channel.empty()
    if (params.gene_bed){
        if (params.gene_bed.endsWith(".gz")) {
            ch_gene_bed = GUNZIP_BED ( params.gene_bed ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    }

    /*
    * Create chromosome sizes file
    */
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES ( ch_fasta ).sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)


    /*
    * Create chromosome sizes file for spike_in
    */
    ch_spikein_chrom_sizes = GET_SPIKEIN_CHROM_SIZES ( ch_spikein_fasta ).sizes

    /*
    * Uncompress Bowtie2 index or generate from scratch if required for both genomes
    */
    ch_bt2_index         = Channel.empty()
    ch_bt2_spikein_index = Channel.empty()
    ch_bt2_versions       = Channel.empty()
    if ("bowtie2" in prepare_tool_indices) {
        if (params.bowtie2) {
            if (params.bowtie2.endsWith(".tar.gz")) {
                ch_bt2_index = UNTAR_BT2_INDEX ( params.bowtie2 ).untar
                ch_versions  = ch_versions.mix(UNTAR_BT2_INDEX.out.versions)
            } else {
                ch_bt2_index = file(params.bowtie2)
            }
        } else {
            ch_bt2_index = BOWTIE2_BUILD ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }

        if (params.spikein_bowtie2) {
            if (params.spikein_bowtie2.endsWith(".tar.gz")) {
                ch_bt2_spikein_index = UNTAR_SPIKEIN_BT2_INDEX ( params.spikein_bowtie2 ).untar
                ch_versions          = ch_versions.mix(UNTAR_SPIKEIN_BT2_INDEX.out.versions)
            } else {
                ch_bt2_spikein_index = file(params.spikein_bowtie2)
            }
        } else {
            ch_bt2_spikein_index = BOWTIE2_SPIKEIN_BUILD ( ch_spikein_fasta ).index
            ch_versions          = ch_versions.mix(BOWTIE2_SPIKEIN_BUILD.out.versions)
        }
    }

    emit:
    fasta                  = ch_fasta                    // path: genome.fasta
    chrom_sizes            = ch_chrom_sizes              // path: genome.sizes
    spikein_chrom_sizes    = ch_spikein_chrom_sizes      // path: genome.sizes
    gtf                    = ch_gtf                      // path: genome.gtf
    bed                    = ch_gene_bed                 // path: genome.bed
    bowtie2_index          = ch_bt2_index                // path: bt2/index/
    bowtie2_spikein_index  = ch_bt2_spikein_index        // path: bt2/index/

    versions               = ch_versions                 // channel: [ versions.yml ]
}
