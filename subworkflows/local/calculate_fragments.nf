/*
 * Calculate bed fragments from bam file
 */

params.samtools_options      = [:]
params.samtools_view_options = [:]
params.bamtobed_options      = [:]
params.awk_options           = [:]
params.cut_options           = [:]

include { SAMTOOLS_VIEW_SORT_STATS } from "./samtools_view_sort_stats"                            addParams( samtools_options: params.samtools_options, samtools_view_options: params.samtools_view_options )
include { BEDTOOLS_BAMTOBED        } from "../../modules/nf-core/software/bedtools/bamtobed/main" addParams( options: params.bamtobed_options )
include { AWK                      } from "../../modules/local/awk"                               addParams( options: params.awk_options      )
include { CUT                      } from "../../modules/local/cut"                               addParams( options: params.cut_options      )

workflow CALCULATE_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    // Filter for mapped reads only
    SAMTOOLS_VIEW_SORT_STATS( bam )

    // Convert to bed file
    BEDTOOLS_BAMTOBED ( SAMTOOLS_VIEW_SORT_STATS.out.bam  )

    // Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    AWK ( BEDTOOLS_BAMTOBED.out.bed  )

    // Only extract the fragment related columns
    CUT ( AWK.out.file )

    emit:
    bed              = CUT.out.file
    mapped_bam       = SAMTOOLS_VIEW_SORT_STATS.out.bam        // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_VIEW_SORT_STATS.out.bai        // channel: [ val(meta), [ bai ] ]
    stats            = SAMTOOLS_VIEW_SORT_STATS.out.stats      // channel: [ val(meta), [ stats ] ]
    flagstat         = SAMTOOLS_VIEW_SORT_STATS.out.flagstat   // channel: [ val(meta), [ flagstat ] ]
    idxstats         = SAMTOOLS_VIEW_SORT_STATS.out.idxstats   // channel: [ val(meta), [ idxstats ] ]
    
    samtools_version = SAMTOOLS_VIEW_SORT_STATS.out.samtools_version //    path: *.version.txt
    bedtools_version = BEDTOOLS_BAMTOBED.out.version                 //    path: *.version.txt
    awk_version      = AWK.out.version                               //    path: *.version.txt
}
