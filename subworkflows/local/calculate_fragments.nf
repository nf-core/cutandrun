/*
 * Calculate bed fragments from bam file
 */

params.samtools_options      = [:]
params.samtools_view_options = [:]
params.samtools_sort_options = [:]
params.bamtobed_options      = [:]
params.awk_options           = [:]
params.cut_options           = [:]

include { SAMTOOLS_VIEW      } from "../../modules/nf-core/modules/samtools/view/main"     addParams( options: params.samtools_view_options )
include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'     addParams( options: params.samtools_sort_options )
include { BEDTOOLS_BAMTOBED  } from "../../modules/nf-core/modules/bedtools/bamtobed/main" addParams( options: params.bamtobed_options      )
include { AWK                } from "../../modules/local/awk"                              addParams( options: params.awk_options           )
include { CUT                } from "../../modules/local/cut"                              addParams( options: params.cut_options           )

workflow CALCULATE_FRAGMENTS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    /*
     * Filter BAM file
     */
    SAMTOOLS_VIEW ( bam )

    /*
    * Sort BAM file
    */
    SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )

    // Convert to bed file
    BEDTOOLS_BAMTOBED ( SAMTOOLS_SORT.out.bam  )

    // Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    AWK ( BEDTOOLS_BAMTOBED.out.bed  )

    // Only extract the fragment related columns
    CUT ( AWK.out.file )

    emit:
    bed              = CUT.out.file          // channel: [ val(meta), [ bed ] ]
    bam              = SAMTOOLS_SORT.out.bam // channel: [ val(meta), [ bam ] ]

    samtools_version = SAMTOOLS_SORT.out.version      //    path: *.version.txt
    bedtools_version = BEDTOOLS_BAMTOBED.out.version  //    path: *.version.txt
    awk_version      = AWK.out.version                //    path: *.version.txt
}
