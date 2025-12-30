/*
 * Pool IgG control BAM files by condition for use with epic2 and SPAN
 */

include { SAMTOOLS_MERGE } from "../../modules/nf-core/samtools/merge/main"
include { SAMTOOLS_INDEX } from "../../modules/nf-core/samtools/index/main"

workflow POOL_CONTROLS {
    take:
    ch_bam_control   // channel: [ val(meta), [ bam ] ]
    pool_by          // value: "condition" or "all"

    main:
    ch_versions = Channel.empty()

    if (pool_by == "all") {
        // Pool all controls together
        ch_bam_control
            .map { meta, bam -> bam }
            .collect()
            .map { bams ->
                def new_meta = [id: "pooled_IgG", condition: "all", is_control: true]
                [new_meta, bams]
            }
            .set { ch_controls_grouped }
    } else {
        // Pool controls by condition (default)
        ch_bam_control
            .map { meta, bam -> [meta.condition, meta, bam] }
            .groupTuple(by: 0)
            .map { condition, metas, bams ->
                def new_meta = [id: "${condition}_pooled_IgG", condition: condition, is_control: true]
                [new_meta, bams.flatten()]
            }
            .set { ch_controls_grouped }
    }

    // Merge BAM files
    SAMTOOLS_MERGE(
        ch_controls_grouped,
        [[],[]],  // fasta (not required for merge)
        [[],[]]   // fai (not required for merge)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Index merged BAM files
    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam      = SAMTOOLS_MERGE.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai      // channel: [ val(meta), [ bai ] ]
    versions = ch_versions                  // channel: [ versions.yml ]
}
