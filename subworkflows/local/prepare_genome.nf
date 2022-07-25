/*
 * Uncompress and prepare reference genome files
*/

include { GUNZIP as GUNZIP_FASTA                     } from '../../modules/nf-core/modules/gunzip/main.nf'
include { GUNZIP as GUNZIP_SPIKEIN_FASTA             } from '../../modules/nf-core/modules/gunzip/main.nf'
include { GUNZIP as GUNZIP_GTF                       } from '../../modules/nf-core/modules/gunzip/main.nf'
include { GUNZIP as GUNZIP_BED                       } from '../../modules/nf-core/modules/gunzip/main.nf'
include { CUSTOM_GETCHROMSIZES as TARGET_CHROMSIZES  } from '../../modules/nf-core/modules/custom/getchromsizes/main.nf'
include { CUSTOM_GETCHROMSIZES as SPIKEIN_CHROMSIZES } from '../../modules/nf-core/modules/custom/getchromsizes/main.nf'
include { UNTAR as UNTAR_INDEX_TARGET                } from '../../modules/nf-core/modules/untar/main.nf'
include { UNTAR as UNTAR_INDEX_SPIKEIN               } from '../../modules/nf-core/modules/untar/main.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_TARGET      } from '../../modules/nf-core/modules/bowtie2/build/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_SPIKEIN     } from '../../modules/nf-core/modules/bowtie2/build/main'
include { TABIX_BGZIPTABIX                           } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { GTF2BED                                    } from '../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX                             } from '../../modules/nf-core/modules/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:
    ch_versions = Channel.empty()

    /*
    * Uncompress genome fasta file if required
    */
    if (params.fasta.endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta    = file(params.fasta)
    }

    /*
    * Uncompress spike-in genome fasta file if required
    */
    if (params.spikein_fasta.endsWith(".gz")) {
        ch_spikein_fasta = GUNZIP_SPIKEIN_FASTA ( [ [:], params.spikein_fasta ] ).gunzip.map { it[1] }
        ch_versions      = ch_versions.mix(GUNZIP_SPIKEIN_FASTA.out.versions)
    } else {
        ch_spikein_fasta = file(params.spikein_fasta)
    }

    /*
    * Uncompress GTF annotation file
    */
    ch_gtf = Channel.empty()
    if (params.gtf.endsWith(".gz")) {
        ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = file(params.gtf)
    }

    ch_gene_bed = Channel.empty()
    if (params.gene_bed){
        /*
        * Uncompress BED annotation file
        */
        if (params.gene_bed.endsWith(".gz")) {
            ch_gene_bed = GUNZIP_BED ( [ [:], params.gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        /*
        * Create GTF bed file if needed
        */
        GTF2BED ( ch_gtf )
        ch_gene_bed = GTF2BED.out.bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    /*
    * Index the bed annotation file
    */
    ch_tabix = ch_gene_bed.map { 
        row -> [ [ id:row.getName() ] , row ] 
    }

    if (params.gene_bed && params.gene_bed.endsWith(".gz")) {
        ch_tabix = ch_tabix.map {
            row ->
                new_id = row[0].id.split("\\.")[0]
                [ [ id: new_id ] , row[1] ] 
        }
    }

    TABIX_BGZIPTABIX (
        ch_tabix
    )
    ch_gene_bed_index = TABIX_BGZIPTABIX.out.gz_tbi
    ch_versions       = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    /*
    * Index genome fasta file
    */
    ch_fasta_index = SAMTOOLS_FAIDX ( ch_fasta.map {row -> [ [:], row ]} ).fai
    ch_versions    = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    /*
    * Create chromosome sizes file
    */
    ch_chrom_sizes = TARGET_CHROMSIZES ( ch_fasta ).sizes
    ch_versions    = ch_versions.mix(TARGET_CHROMSIZES.out.versions)

    /*
    * Create chromosome sizes file for spike_in
    */
    ch_spikein_chrom_sizes = SPIKEIN_CHROMSIZES ( ch_spikein_fasta ).sizes

    /*
    * Uncompress Bowtie2 index or generate from scratch if required for both genomes
    */
    ch_bt2_index         = Channel.empty()
    ch_bt2_spikein_index = Channel.empty()
    ch_bt2_versions      = Channel.empty()
    if ("bowtie2" in prepare_tool_indices) {
        if (params.bowtie2) {
            if (params.bowtie2.endsWith(".tar.gz")) {
                ch_bt2_index = UNTAR_INDEX_TARGET ( [ [], params.bowtie2 ] ).untar.map{ row -> [ row[1] ] }
                ch_versions  = ch_versions.mix(UNTAR_INDEX_TARGET.out.versions)
            } else {
                ch_bt2_index = file(params.bowtie2)
            }
        } else {
            ch_bt2_index = BOWTIE2_BUILD_TARGET ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BOWTIE2_BUILD_TARGET.out.versions)
        }

        if (params.spikein_bowtie2) {
            if (params.spikein_bowtie2.endsWith(".tar.gz")) {
                ch_bt2_spikein_index = UNTAR_INDEX_SPIKEIN ( [ [], params.spikein_bowtie2 ] ).untar.map{ row -> [ row[1] ] }
                ch_versions          = ch_versions.mix(UNTAR_INDEX_SPIKEIN.out.versions)
            } else {
                ch_bt2_spikein_index = file(params.spikein_bowtie2)
            }
        } else {
            ch_bt2_spikein_index = BOWTIE2_BUILD_SPIKEIN ( ch_spikein_fasta ).index
            ch_versions          = ch_versions.mix(BOWTIE2_BUILD_SPIKEIN.out.versions)
        }
    }

    emit:
    fasta                  = ch_fasta                    // path: genome.fasta
    fasta_index            = ch_fasta_index              // path: genome.fai
    chrom_sizes            = ch_chrom_sizes              // path: genome.sizes
    spikein_chrom_sizes    = ch_spikein_chrom_sizes      // path: genome.sizes
    gtf                    = ch_gtf                      // path: genome.gtf
    bed                    = ch_gene_bed                 // path: genome.bed
    bed_index              = ch_gene_bed_index           // path: genome.bed_index

    bowtie2_index          = ch_bt2_index                // path: bt2/index/
    bowtie2_spikein_index  = ch_bt2_spikein_index        // path: bt2/index/

    versions               = ch_versions                 // channel: [ versions.yml ]
}
