/*
 * Perform full suite of deep tools analysis on bam files
*/

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../modules/nf-core/deeptools/multibamsummary/main'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA         } from '../../modules/nf-core/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTFINGERPRINT } from '../../modules/nf-core/deeptools/plotfingerprint/main'

workflow DEEPTOOLS_QC {
    take:
    bam         // channel: [ val(meta), [ bam ] ]
    bai         // channel: [ val(meta), [ bai ] ]
    corr_method // val

    main:
    ch_versions = Channel.empty()

    /*
    * CHANNEL: Combine bam and bai files on id
    */
    bam
    .join( bai )
    .set { ch_bam_bai }
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    // ch_bam_bai | view

    /*
    * CHANNEL: Get list of sample ids
    */
    ch_bam_bai.map { row -> [row[0].id] }
    .collect()
    .map { row -> [row] }
    .set { ch_ids }
    //ch_ids | view

    /*
    * CHANNEL: Combine bam and bai files into one list
    * if we only have one file then cancel correlation and PCA
    */
    ch_bam_bai.map { row -> [row[1]] }
    .collect()
    .map { row -> [row] }
    .combine( ch_bam_bai.map { row -> [row[2]] }.collect().map { row -> [row] } )
    .combine( ch_ids )
    .map { row -> [[id: 'all_target_bams'], row[0], row[1], row[2], row[1].size()] }
    .filter { row -> row[4] > 1 }
    .map { row -> [row[0], row[1], row[2], row[3]] }
    .set { ch_bam_bai_all }
    //ch_bam_bai_all | view

    /*
    * MODULE: Summarise bams into bins
    */
    DEEPTOOLS_MULTIBAMSUMMARY (
        ch_bam_bai_all
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_MULTIBAMSUMMARY.out.versions)
    //DEEPTOOLS_MULTIBAMSUMMARY.out.matrix | view

    /*
    * MODULE: Plot correlation matrix
    */
    DEEPTOOLS_PLOTCORRELATION (
        DEEPTOOLS_MULTIBAMSUMMARY.out.matrix,
        corr_method,
        "heatmap"
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCORRELATION.out.versions)
    //DEEPTOOLS_MULTIBAMSUMMARY.out.matrix | view

    /*
    * MODULE: Plot PCA's
    */
    DEEPTOOLS_PLOTPCA (
        DEEPTOOLS_MULTIBAMSUMMARY.out.matrix
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPCA.out.versions)
    //DEEPTOOLS_PLOTPCA.out.matrix | view

    /*
    * MODULE: Plot Fingerprint
    */
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions)
    //DEEPTOOLS_PLOTFINGERPRINT.out.matrix | view

    emit:
    correlation_matrix = DEEPTOOLS_PLOTCORRELATION.out.matrix
    pca_data           = DEEPTOOLS_PLOTPCA.out.tab
    fingerprint_matrix = DEEPTOOLS_PLOTFINGERPRINT.out.matrix

    versions               = ch_versions                 // channel: [ versions.yml ]
}
