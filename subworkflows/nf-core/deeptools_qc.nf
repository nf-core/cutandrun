/*
 * Perform full suite of deep tools analysis on bam files
*/

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../modules/local/deeptools/multibamsummary/main'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/local/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA         } from '../../modules/local/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTFINGERPRINT } from '../../modules/nf-core/modules/deeptools/plotfingerprint/main'

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
    * CHANNEL: Get list of sample ids
    */
    ch_bam_target.map { row -> [row[0].id] }
    .collect()
    .map { row -> [row] }
    .set { ch_ids }
    //ch_ids | view

    /*
    * CHANNEL: Combine bam and bai files into one list
    * if we only have one file then cancel correlation and PCA
    */
    ch_bam_target.map { row -> [row[1]] }
    .collect()
    .map { row -> [row] } 
    .combine( ch_bai_target.map { row -> [row[1]] }.collect().map { row -> [row] } )
    .combine( ch_ids )
    .map { row -> [[id: 'all_target_bams'], row[0], row[1], row[2], row[1].size()] }
    .filter { row -> row[4] > 1 }
    .map { row -> [row[0], row[1], row[2], row[3]] }
    .set { ch_bam_bai_all }
    ch_bam_bai_all | view

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
        DEEPTOOLS_MULTIBAMSUMMARY.out.matrix
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
