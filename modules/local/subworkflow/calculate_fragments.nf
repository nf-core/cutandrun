/*
 * Calculate bed fragments from bam file
 */

params.samtools_options         = [:]
params.samtools_view_options    = [:]
params.bamtobed_options    = [:]
params.awk_options    = [:]
params.cut_options    = [:]

include { SAMTOOLS_VIEW_SORT_STATS } from './samtools_view_sort_stats' addParams( samtools_options: params.samtools_options , samtools_view_options: params.samtools_view_options)
include { BEDTOOLS_BAMTOBED } from '../process/bedtools_bamtobed' addParams( options: params.bamtobed_options)
include { AWK } from '../process/awk' addParams( options: params.awk_options)
include { CUT } from '../process/cut' addParams( options: params.cut_options)

workflow CALCULATE_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    // Filer for mapped reads only
    SAMTOOLS_VIEW_SORT_STATS( bam )

    // Convert to bed file
    BEDTOOLS_BAMTOBED ( SAMTOOLS_VIEW_SORT_STATS.out.bam  )

    // Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    AWK ( BEDTOOLS_BAMTOBED.out.bed  )

    // Only extract the fragment related columns
    CUT ( AWK.out.file )

    emit:
    bed              = CUT.out.file
    mapped_bam       = SAMTOOLS_VIEW_SORT_STATS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_VIEW_SORT_STATS.out.bai            // channel: [ val(meta), [ bai ] ]
    stats            = SAMTOOLS_VIEW_SORT_STATS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat         = SAMTOOLS_VIEW_SORT_STATS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats         = SAMTOOLS_VIEW_SORT_STATS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    
    samtools_version = SAMTOOLS_VIEW_SORT_STATS.out.samtools_version        //    path: *.version.txt
    bedtools_version = BEDTOOLS_BAMTOBED.out.version        //    path: *.version.txt
    awk_version      = AWK.out.version        //    path: *.version.txt
}






// params.align_options            = [:]
// params.spikein_align_options    = [:]
// params.samtools_options         = [:]
// params.samtools_spikein_options = [:]

// include { BOWTIE2_ALIGN                                  } from '../software/bowtie2/align/main'              addParams( options: params.align_options, save_unaligned: params.save_unaligned )
// include { BOWTIE2_ALIGN as BOWTIE2_SPIKEIN_ALIGN         } from '../software/bowtie2/align/main'              addParams( options: params.spikein_align_options, save_unaligned: false )
// include { BAM_SORT_SAMTOOLS                              } from '../../nf-core/subworkflow/bam_sort_samtools' addParams( options: params.samtools_options )
// include { BAM_SORT_SAMTOOLS as BAM_SORT_SAMTOOLS_SPIKEIN } from '../../nf-core/subworkflow/bam_sort_samtools' addParams( options: params.samtools_spikein_options )

// workflow ALIGN_BOWTIE2 {
//     take:
//     reads         // channel: [ val(meta), [ reads ] ]
//     index         // channel: /path/to/bowtie2/target/index/
//     spikein_index // channel: /path/to/bowtie2/spikein/index/
    
//     main:
//     /*
//      * Map reads with BOWTIE2 to target genome
//      */
//     BOWTIE2_ALIGN ( reads, index )

//     /*
//      * Map reads with BOWTIE2 to spike-in genome
//      */
//     BOWTIE2_SPIKEIN_ALIGN ( reads, spikein_index )

//     /*
//      * Sort, index BAM file and run samtools stats, flagstat and idxstats
//      */
//     BAM_SORT_SAMTOOLS         ( BOWTIE2_ALIGN.out.bam         )
//     BAM_SORT_SAMTOOLS_SPIKEIN ( BOWTIE2_SPIKEIN_ALIGN.out.bam )

//     emit:
//     bowtie2_version     = BOWTIE2_ALIGN.out.version              // path: *.version.txt
//     samtools_version    = BAM_SORT_SAMTOOLS.out.version           //    path: *.version.txt
    
//     orig_bam            = BOWTIE2_ALIGN.out.bam                  // channel: [ val(meta), bam            ]
//     orig_spikein_bam    = BOWTIE2_SPIKEIN_ALIGN.out.bam          // channel: [ val(meta), bam            ]
    
//     bowtie2_log         = BOWTIE2_ALIGN.out.log                  // channel: [ val(meta), log_final      ]
//     bowtie2_spikein_log = BOWTIE2_SPIKEIN_ALIGN.out.log          // channel: [ val(meta), log_final      ]
    
//     bam                 = BAM_SORT_SAMTOOLS.out.bam              // channel: [ val(meta), [ bam ] ]
//     bai                 = BAM_SORT_SAMTOOLS.out.bai              // channel: [ val(meta), [ bai ] ]
//     stats               = BAM_SORT_SAMTOOLS.out.stats            // channel: [ val(meta), [ stats ] ]
//     flagstat            = BAM_SORT_SAMTOOLS.out.flagstat         // channel: [ val(meta), [ flagstat ] ]
//     idxstats            = BAM_SORT_SAMTOOLS.out.idxstats         // channel: [ val(meta), [ idxstats ] ]

//     spikein_bam         = BAM_SORT_SAMTOOLS_SPIKEIN.out.bam      // channel: [ val(meta), [ bam ] ]
//     spikein_bai         = BAM_SORT_SAMTOOLS_SPIKEIN.out.bai      // channel: [ val(meta), [ bai ] ]
//     spikein_stats       = BAM_SORT_SAMTOOLS_SPIKEIN.out.stats    // channel: [ val(meta), [ stats ] ]
//     spikein_flagstat    = BAM_SORT_SAMTOOLS_SPIKEIN.out.flagstat // channel: [ val(meta), [ flagstat ] ]
//     spikein_idxstats    = BAM_SORT_SAMTOOLS_SPIKEIN.out.idxstats // channel: [ val(meta), [ idxstats ] ]
// }
