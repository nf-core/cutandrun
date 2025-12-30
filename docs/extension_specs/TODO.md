# 1. Add "condition" to samplesheet (ie. Control | Treatment}
# 2. Add group-specific handling for:
    - Normalization
    - Peak-Calling
# 3. Add support for additional peakcallers:
    - GoPeaks
    - Epic2
    - SPAN / OMNIPEAKS
    - MACS2 narrow + broad
    - GoPeaks Narrow + broad
# 4. Add differential peak calling sub workflow
    - diffbind (generate samplesheets from called peaks)
    - ChIPBinner (standalone)
    - SPAN (standalone)
# 5. Add --compare_norm_methods mode
# 6. Add --compare_peak_callers mode
# 7. Add differential peak / enrichment characterization
    - ChromHMM
    - cCRE
    - {other features?}
    - Create LOLA/GSEA sets from group x cond
