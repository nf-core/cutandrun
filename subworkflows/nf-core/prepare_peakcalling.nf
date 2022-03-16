/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { BEDTOOLS_GENOMECOV    } from "../../modules/nf-core/modules/bedtools/genomecov/main"
include { DEEPTOOLS_BAMCOVERAGE } from "../../modules/local/modules/deeptools/bamcoverage/main"
include { BEDTOOLS_SORT         } from "../../modules/nf-core/modules/bedtools/sort/main"
include { UCSC_BEDCLIP          } from "../../modules/nf-core/modules/ucsc/bedclip/main"
include { UCSC_BEDGRAPHTOBIGWIG } from "../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main"

workflow PREPARE_PEAKCALLING {
    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_chrom_sizes // channel: [ sizes ]
    ch_dummy_file  // channel: [ dummy ]
    norm_mode      // value: ["Spikein", "RPKM", "CPM", "BPM", "RPGC", "None" ]

    main:
    ch_versions = Channel.empty()
    ch_bam_out  = Channel.empty()
    ch_bedgraph = Channel.empty()

    if (norm_mode == "Spikein") {
        /*
        * CHANNEL: Calculate scale factor for each sample based on a constant devided by the number
        *          of reads aligned to the spike-in genome.
        */
        ch_bam
            .map { row ->
                def denominator = row[0].find{ it.key == "bt2_total_aligned_spikein" }?.value.toInteger()
                [ row[0].id, params.normalisation_c / (denominator != 0 ? denominator : 1) ]
            }
        .set { ch_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [id, scale_factor]
        //ch_scale_factor | view
    }
    else if (norm_mode == "None") {
        /*
        * CHANNEL: Assign scale factor of 1
        */
        ch_bam
            .map { row ->
                [ row[0].id, 1 ]
            }
        .set { ch_scale_factor }
    }

    if (norm_mode == "Spikein" || norm_mode == "None") {
        /*
        * CHANNEL: Create a channel with the scale factor as a seperate value
        */
        ch_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_scale_factor )
            .map { row -> row[1..(row.size() - 1)] }
            .map { row ->
                row[0].put("scale_factor", row[2])
                [ row[0], row[1], row[2] ] }
        .set { ch_bam_scale }
        //EXAMPLE CHANNEL STRUCT: [[META + scale_factor:10000], BAM, SCALE_FACTOR]
        //ch_bam_scale | view

        /*
        * MODULE: Convert bam files to bedgraph
        */
        BEDTOOLS_GENOMECOV (
            ch_bam_scale,
            ch_dummy_file,
            "bedGraph"
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
        ch_bedgraph = BEDTOOLS_GENOMECOV.out.genomecov
        //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
        //BEDTOOLS_GENOMECOV.out.genomecov | view

        /*
        * CHANNEL: Add the scale factor values to the main meta-data stream
        */
        ch_bam_scale
            .map { row -> [ row[0], row[1] ] }
        .set { ch_bam_out }
        //EXAMPLE CHANNEL STRUCT: [[META], BAM]
        //ch_samtools_bam | view
    } else {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bam_bai | view

        /*
        * CHANNEL: Split files based on igg or not
        */
        ch_bam_bai.branch { it ->
            target: it[0].group != "igg"
            control: it[0].group == "igg"
        }
        .set { ch_bam_bai_split }

        /*
        * CHANNEL: Assign scale factor of 1 to target files
        */
        ch_bam_bai_split.target
            .map { row ->
                [ row[0], row[1], row[2], 1 ]
            }
        .set { ch_bam_bai_split_target }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_target | view

        /*
        * CHANNEL: Assign igg scale factor to target files
        */
        ch_bam_bai_split.control
            .map { row ->
                [ row[0], row[1], row[2], params.igg_scale_factor ]
            }
        .set { ch_bam_bai_split_igg }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_igg | view

        /*
        * CHANNEL: Mix the split channels back up
        */
        ch_bam_bai_split_target
            .mix(ch_bam_bai_split_igg)
        .set { ch_bam_bai_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_scale_factor | view

        /*
        * MODULE: Convert bam files to bedgraph and normalise
        */
        DEEPTOOLS_BAMCOVERAGE (
            ch_bam_bai_scale_factor
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)
        ch_bedgraph = DEEPTOOLS_BAMCOVERAGE.out.bedgraph
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bedgraph | view

        // Dont assign any new meta data
        ch_bam_out = ch_bam
    }

    /*
    * MODULE: Sort bedgraph
    */
    BEDTOOLS_SORT (
        ch_bedgraph,
        "bedGraph"
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    /*
    * MODULE: Clip off bedgraphs so none overlap beyond chromosome edge
    */
    UCSC_BEDCLIP (
        BEDTOOLS_SORT.out.sorted,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions)
    //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
    //UCSC_BEDCLIP.out.bedgraph | view

    /*
    * MODULE: Convert bedgraph to bigwig
    */
    UCSC_BEDGRAPHTOBIGWIG (
        UCSC_BEDCLIP.out.bedgraph,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], BIGWIG]
    //UCSC_BEDGRAPHTOBIGWIG.out.bigwig | view

    emit:
    bam      = ch_bam_out                       // channel: [ val(meta), [ bam ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
