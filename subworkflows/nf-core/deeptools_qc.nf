/*
 * Perform full suite of deep tools analysis on bam files
*/

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../modules/local/deeptools/multibamsummary/main.nf'

workflow DEEPTOOLS_QC {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    bai // channel: [ val(meta), [ bai ] ]

    main:
    ch_versions = Channel.empty()

    /*
    * CHANNEL: Filter bams for target only
    */
    bam.filter { it -> it[0].is_control == false }
    .set { ch_bam_target }
    //ch_bam_target | view

    /*
    * CHANNEL: Filter bais for target only
    */
    bai.filter { it -> it[0].is_control == false }
    .set { ch_bai_target }
    //ch_bai_target | view

    /*
    * CHANNEL: Combine bam and bai files on id
    */
    ch_bam_target.map { row -> [row[0].id, row ].flatten()}
    .join ( ch_bai_target.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [row[1], row[2], row[4]] }
    .set { ch_bam_bai }
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bam_bai | view

    /*
    * CHANNEL: Combine bam and bai files into one list
    */
    ch_bam_target.map { row -> [row[1]] }
    .collect()
    .map { row -> [row] } 
    .combine( ch_bai_target.map { row -> [row[1]] }.collect().map { row -> [row] } )
    .map { row -> [[:], row[0], row[1]] } 
    .set { ch_bam_bai_all }
    //ch_bam_bai_all | view

    /*
    * MODULE: Summarise bams into bins
    */
    DEEPTOOLS_MULTIBAMSUMMARY (
        ch_bam_bai_all
    )
    //DEEPTOOLS_MULTIBAMSUMMARY.out.matrix | view


    emit:
    // fasta                  = ch_fasta                    // path: genome.fasta
    // fasta_index            = ch_fasta_index              // path: genome.fai
    // chrom_sizes            = ch_chrom_sizes              // path: genome.sizes
    // spikein_chrom_sizes    = ch_spikein_chrom_sizes      // path: genome.sizes
    // gtf                    = ch_gtf                      // path: genome.gtf
    // bed                    = ch_gene_bed                 // path: genome.bed
    // bed_index              = ch_gene_bed_index           // path: genome.bed_index
    // allowed_regions        = ch_genome_include_regions   // path: genome.regions

    // bowtie2_index          = ch_bt2_index                // path: bt2/index/
    // bowtie2_spikein_index  = ch_bt2_spikein_index        // path: bt2/index/

    versions               = ch_versions                 // channel: [ versions.yml ]
}
