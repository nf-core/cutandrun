/*
 * First use custom .py script to find unique linear amplification alignments, samtools filter unique alignments, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BEDTOOLS_BAMTOBED       } from "../../modules/nf-core/bedtools/bamtobed/main"
include { FIND_UNIQUE_READS       } from '../../modules/local/python/find_unique_reads'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_SORT           } from "../../modules/nf-core/samtools/sort/main.nf"
include { SAMTOOLS_VIEW           } from "../../modules/nf-core/samtools/view/main.nf"

workflow DEDUPLICATE_LINEAR {
    take:
    bam            // channel: [ val(meta), [ bam ] ]
    bai            // channel: [ val(meta), [ bai ] ]
    fasta          // channel: [ val(meta), fasta ]
    fai            // channel: [ val(meta), fai ]
    process_target // boolean
    mqc_header     // path

    main:
    /*
    * Find unique linear amplification alignments
    */
    ch_bam           = Channel.empty()
    ch_metrics       = Channel.empty()
    ch_versions      = Channel.empty()
    ch_linear_duplicates = Channel.empty()

    if( process_target ) {

        SAMTOOLS_SORT (
            bam
        )
        ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

        // Convert .bam files to bed to find unique Tn5-ME-A insertion sites
        BEDTOOLS_BAMTOBED (
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix( BEDTOOLS_BAMTOBED.out.versions )

        // Use custom .py script to find names of unique alignments
        FIND_UNIQUE_READS (
            BEDTOOLS_BAMTOBED.out.bed,
            mqc_header
        )
        ch_linear_duplicates = FIND_UNIQUE_READS.out.txt
        ch_metrics           = FIND_UNIQUE_READS.out.metrics
        ch_versions          = ch_versions.mix( FIND_UNIQUE_READS.out.versions )

        // Subset original .bam file to contain only unique alignments
        bam
        .join( bai )
        .join( ch_linear_duplicates )
        .set { ch_bam_bai_txt }

        ch_bam_bai_txt
        .map { row -> row[0,1,2]}
        .set { ch_bam_bai }

        ch_bam_bai_txt
        .map { row -> row[3]}
        .set { ch_txt}

        SAMTOOLS_VIEW (
            ch_bam_bai,
            fasta,
            ch_txt
        )
        ch_versions = ch_versions.mix( SAMTOOLS_VIEW.out.versions )
        ch_bam      = SAMTOOLS_VIEW.out.bam

    }
    else { // Split out control files and run only on these

        bam.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_split }
        //ch_split.target | view
        //ch_split.control | view

        SAMTOOLS_SORT (
            ch_split.control
        )
        ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

        // Run .bam to bed on control files only and find unique alignments in control files
        BEDTOOLS_BAMTOBED (
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix( BEDTOOLS_BAMTOBED.out.versions )

        // Use custom .py script to find names of unique alignments in control files only
        FIND_UNIQUE_READS (
            BEDTOOLS_BAMTOBED.out.bed,
            mqc_header
        )
        ch_linear_duplicates = FIND_UNIQUE_READS.out.txt
        ch_metrics           = FIND_UNIQUE_READS.out.metrics
        ch_versions          = ch_versions.mix( FIND_UNIQUE_READS.out.versions )

        // Subset original .bam file to contain only unique alignments
        ch_split.control
        .join( bai )
        .join( ch_linear_duplicates )
        .set { ch_bam_bai_txt }

        ch_bam_bai_txt
        .map { row -> row[0,1,2]}
        .set { ch_bam_bai }

        ch_bam_bai_txt
        .map { row -> row[3]}
        .set { ch_txt}

        // Subset original .bam file to contain only unique alignments
        SAMTOOLS_VIEW (
            ch_bam_bai,
            fasta,
            ch_txt
        )
        ch_versions = ch_versions.mix( SAMTOOLS_VIEW.out.versions )
        ch_bam = SAMTOOLS_VIEW.out.bam

        // Prevents issues with resume with the branch elements coming in the wrong order
        ch_sorted_targets = ch_split.target
            .toSortedList { row -> row[0].id }
            .flatMap()

        ch_sorted_controls = ch_bam
            .toSortedList { row -> row[0].id }
            .flatMap()

        // Return the filtered control bam file concatenated with the original target bam
        ch_bam = ch_sorted_targets.concat( ch_sorted_controls )
    }

    /*
    * WORKFLOW: Re sort and index all the bam files + calculate stats
    */
    BAM_SORT_STATS_SAMTOOLS (
        ch_bam,
        fasta
    )
    ch_versions = ch_versions.mix( BAM_SORT_STATS_SAMTOOLS.out.versions )

    emit:
    bam                = BAM_SORT_STATS_SAMTOOLS.out.bam           // channel: [ val(meta), [ bam ] ]
    bai                = BAM_SORT_STATS_SAMTOOLS.out.bai           // channel: [ val(meta), [ bai ] ]
    stats              = BAM_SORT_STATS_SAMTOOLS.out.stats         // channel: [ val(meta), [ stats ] ]
    flagstat           = BAM_SORT_STATS_SAMTOOLS.out.flagstat      // channel: [ val(meta), [ flagstat ] ]
    idxstats           = BAM_SORT_STATS_SAMTOOLS.out.idxstats      // channel: [ val(meta), [ idxstats ] ]
    metrics            = ch_metrics                                // channel: [ metrics.txt  ]
    linear_metrics_mqc = FIND_UNIQUE_READS.out.linear_metrics_mqc  // channel: [ mqc.tsv ]
    versions           = ch_versions                               // channel: [ versions.yml ]
}
