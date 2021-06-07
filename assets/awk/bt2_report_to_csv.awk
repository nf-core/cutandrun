#!/bin/awk -f

# This awk script parses a bowtie2 alignment summary to a csv file for easier downstream processing

BEGIN {
    cols="bt2_total_reads:bt2_align1:bt2_align_gt1:bt2_non_aligned:bt2_total_aligned"
    FS="\n"
    RS="\n\n"
    col_count=split(cols, col_arr, ":");
    for(i=1; i<=col_count; i++) printf col_arr[i] ((i==col_count) ? "\n" : ",");
}
{
    for (i=1; i<=NF; i++) {
        if(index($i,"reads; of these:") != 0) {
            split($i, line_split, " ");
            data["bt2_total_reads"]=line_split[1];
        }
        if(index($i,"aligned concordantly exactly 1 time") != 0) {
            split($i, line_split, " ");
            data["bt2_align1"]=line_split[1];
        }
        if(index($i,"aligned concordantly >1 times") != 0) {
            split($i, line_split, " ");
            data["bt2_align_gt1"]=line_split[1];
        }
        if(index($i,"aligned concordantly 0 times") != 0) {
            split($i, line_split, " ");
            data["bt2_non_aligned"]=line_split[1];
        }
    }
}
END {
    data["bt2_total_aligned"] = data["bt2_align1"] + data["bt2_align_gt1"]
    for (i=1; i<=col_count; i++) printf data[col_arr[i]] ((i==col_count) ? "\n" : ",");
}
