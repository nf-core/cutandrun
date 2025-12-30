# 1. Add "condition" to samplesheet (ie. Control | Treatment}
# 2. Add support for additional peakcallers:
    - GoPeaks
    - Epic2
    - SPAN / OMNIPEAKS
    - MACS2 narrow + broad
    - GoPeaks Narrow + broad
# 3. Add differential peak calling sub workflow
    - diffbind (generate samplesheets from called peaks)
    - ChIPBinner (standalone)
    - SPAN (standalone)
# 4. Add --compare_norm_methods mode
# 5. Add --compare_peak_callers mode
# 6. Add differential peak / enrichment characterization
    - ChromHMM
    - cCRE
    - {other features?}
    - Create LOLA/GSEA sets from group x cond
