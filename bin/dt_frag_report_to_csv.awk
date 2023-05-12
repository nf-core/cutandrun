#!/bin/awk -f

BEGIN {
    cols="dt_frag_sampled:dt_frag_mean_len:dt_frag_min_len:dt_frag_max_len"
    FS="\t"
    RS="\n"
    col_count=split(cols, col_arr, ":");
    for(i=1; i<=col_count; i++) printf col_arr[i] ((i==col_count) ? "\n" : ",");
}
{
    for (i=1; i<=NF; i++) {
        if(index($i,"Frag. Sampled") !=0) {
            data["dt_frag_sampled"]=i;
        }
        if(index($i,"Frag. Len. Mean") !=0) {
            data["dt_frag_mean_len"]=i;
        }
        if(index($i,"Frag. Len. Min") !=0) {
            data["dt_frag_min_len"]=i;
        }
        if(index($i,"Frag. Len. Max") !=0) {
            data["dt_frag_max_len"]=i;
        }
    }
}
END {
    for (j=1; j<=col_count; j++) printf $(data[col_arr[j]]) ((j==col_count) ? "\n" : ",");
}
