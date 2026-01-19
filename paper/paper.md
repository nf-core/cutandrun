---
title: 'nf-core/cutandrun: A Nextflow pipeline for the analysis of CUT&RUN, CUT&Tag and TIP-seq datasets'
tags:
  - Nextflow
  - nf-core
  - CUT&RUN
  - CUT&Tag
  - TIP-seq
  - single-cell
  - epigenomics
authors:
  - name: Tamara L. Hodgetts
    orcid: 0000-0002-3991-1683
    corresponding: true
    affiliation: 1 
  - name: Charlotte West
    orcid: 0000-0003-2094-3324
    affiliation: 2
  - name: Nicholas M. Luscombe
    orcid: 0000-0001-5293-4778
    affiliation: 3
  - name: Jernej Ule
    orcid: 0000-0002-2452-4277
    affiliation: 1
  - name: James Briscoe
    orcid: 0000-0002-1020-5240
    affiliation: 1
  - name: Chris Cheshire
    orcid: 0000-0002-5806-2883
    affiliation: 1

affiliations:
 - name: The Francis Crick Institute, London, UK
   index: 1
 - name: European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI), Wellcome Genome Campus, Hinxton, UK
   index: 2
 - name: Okinawa Institute of Science and Technology, Okinawa, Japan
   index: 3
date: 7 March 2025
bibliography: paper.bib
---

# Summary

Mapping transcription factor binding events and histone modifications on a genome-wide scale is key to understanding the regulatory dynamics controlling gene expression. Recently, CUT&RUN, CUT&Tag and TIP-seq protocols provide new avenues to map chromatin-associated proteins, offering advantages compared to ChIP-seq in terms of sensitivity, resolution and sample consumption. Despite differences in biochemistry, the core computational analysis of CUT&RUN, CUT&Tag and TIP-seq data exhibit notable similarities, representing an opportunity for a unified bioinformatics solution. 

Here, we present nf-core/cutandrun, a best-practice bioinformatics analysis workflow for CUT&RUN, CUT&Tag and TIP-seq data. In contrast to existing alternative pipelines, nf-core/cutandrun was developed based on the nf-core community framework to provide enhanced reproducibility, flexibility, scalability, portability and robustness. nf-core/cutandrun additionally enables the specification of spike-in genomes with different normalisation options, supports multiple peak callers, and provides extensive quality control metrics reporting. The pipeline supports a wide range of execution environments and operating systems, and is usable with minimal technical knowledge. We encourage user engagement with the wider nf-core community through the official nf-core GitHub repository and Slack channels, and aim to incorporate feedback in subsequent releases to ensure that the pipeline remains dynamic and up-to-date.


# Statement of need

Characterising the genome-wide distribution of transcription factor binding events and histone modifications is key to elucidating the mechanics and dynamics of gene expression regulation. In recent years, the CUT&RUN and CUT&Tag protocols have emerged as alternatives to chromatin immunoprecipitation (ChIP), providing enhanced signal-to-noise ratios, as well as reduced costs and input sample requirements [@Skene:2017; @Kaya-Okur:2019; see for a review @Furey:2012; @Park:2009]. The CUT&RUN and CUT&Tag protocols first expose permeabilised, unfixed cells to primary antibodies, which bind a particular DNA-associated transcription factor or histone modification of interest. Secondary antibodies are added for signal amplification, and to facilitate the binding of either Protein A-MNase (CUT&RUN) or Protein A-Tn5 transposase (CUT&Tag) fusion constructs. Enzymatic activation subsequently enables local DNA cleavage to generate protein-associated chromatin fragments, which is followed by DNA extraction, library preparation, amplification and sequencing. Most recently, the TIP-seq protocol was additionally developed, which combines CUT&Tag with *in vitro* transcription to provide a further 10-fold increase in sensitivity [@Bartlett:2021]. 

Here, we present nf-core/cutandrun, a best-practice bioinformatics pipeline for CUT&RUN, CUT&Tag and TIP-seq datasets. To meet the modern requirements for bioinformatics analysis, which we define to be reproducibility, flexibility, scalability, portability and robust reporting, we developed the pipeline in line with the nf-core community framework [@Ewels:2020; @Tommaso:2017]. nf-core/cutandrun provides several advantages relative to CUT&RUNTools and ePeak, which are the only current alternative pipelines for the end-to-end analysis of CUT&RUN and CUT&Tag data [@Zhu:2019; @Yu:2022; @Daunesse:2022]. First, as a result of the protein production process in original CUT&RUN and CUT&Tag protocols, there is residual *E. coli* DNA that can be used as a natural spike-in calibrator. This *E. coli* DNA is increasingly purified out of commercial kits, with spike-in DNA being manually added during the protocol. Manual addition of spike-in DNA can introduce considerable experimental error, leading to low alignment rates that are insufficient for normalisation. Consequently, we included additional read count-based normalisation modes and support for user-defined spike-in genomes. The pipeline additionally provides three distinct peak calling methods, including SEACR, MACS2 narrow, and MACS2 broad peak modes. We selected SEACR as the default peak caller, as it has been demonstrated to provide greater specificity and efficiency of peak calling on sparse data [@Meers:2019]. Finally, nf-core/cutandrun provides extensive quality control (QC) metrics reporting to maximise confidence prior to downstream analysis, and uniquely supports the analysis of TIP-seq data. 


# Workflow overview

nf-core/cutandrun offers a robust data analysis workflow that is designed for Illumina-based sequencing data obtained from CUT&RUN, CUT&Tag and TIP-seq experiments (Figure 1).

![A schematic representation of the nf-core/cutandrun pipeline.
Icons were adapted from [nf-core/sarek](https://nf-core/sarek); under a [CC-BY 4.0 licence](https://creativecommons.org/licenses/by/4.0/).](figure1.png)

The workflow initialises by performing checks on input files and determining the experimental design in terms of sample grouping, potential IgG controls, and the required normalisation mode. During pre-processing, the raw FASTQ files are merged, considering multi-lane sequencing or technical replicates as defined in the sample sheet. The resultant, merged FASTQ files are initially screened using FastQC to collect QC metrics, before undergoing trimming by TrimGalore to remove sequencing adaptors from reads [@Andrews:2010; @Krueger:2012]. FastQC is re-executed post-trimming, consistent with the workflow’s emphasis on transparency and interpretability. The reads are next aligned to both the target reference genome and the spike-in genome using Bowtie2 [@Langmead:2012]. Alignment statistics, fragment distributions and binned sample correlation metrics are calculated for each sample using Samtools, deepTools, and custom Python scripts, aiding users in identifying potential issues with the alignment process [@Li:2009; @Ramirez:2016]. Reads which aligned to the target genome are subjected to q-score filtering using Samtools, and deduplication using Picard (IgG controls only by default) [@Li:2009; @Picard:2019]. For TIP-seq data, linear amplification duplicates are additionally removed with custom Python scripts. Scaling factors for BAM-to-BedGraph conversion are calculated on a per-sample basis, based on the number of reads that aligned to the spike-in genome. BedGraph files are also converted to a BigWig format to auto-generate IGV sessions, and generate heatmaps using deepTools [@Robinson:2011; @Ramirez:2016]. By default, peak calling is performed using SEACR, which was specifically developed for protocols with low background noise [@Meers:2019]. However, the pipeline also supports peak calling using MACS2 either instead of or in addition to SEACR. At this stage, any IgG controls are split out and used to normalise peak calling against background levels of non-specific antibody binding. Peaks are merged into consensus peak sets based on sample groupings, and intra-group peak reproducibility is assessed using both the number of reproducible peaks and FRiP scores (Fraction of Reads in Peaks). Finally, comprehensive summaries of the workflow execution, parameter settings, software versions, and QC metrics are aggregated into a single MultiQC report [@Ewels:2016].


# Pipeline installation and operation 

nf-core/cutandrun offers minimal execution requirements, requiring only the installation of Nextflow (version $\geq$ 21.10.3) on the host system and a containerisation platform such as Docker [see for a review @Rad:2017] or Singularity [@Kurtzer:2017]. The workflow is automatically downloaded by Nextflow when the keywords 'nf-core/cutandrun' are included within a run statement, rendering manual installation unnecessary. Pipeline installation and operation are easily verified by running pre-defined test routines with public data hosted in the cloud, allowing users to stress test the host system before data analysis. nf-core/cutandrun supports a wide range of execution environments including local, cluster and cloud instances that run on macOS, Linux or Windows operating systems. A rich set of parameters is provided to support workflow customisation and flow-switching, and all results are provided in a single, intuitively-structured folder to facilitate troubleshooting. Readers are referred to the official [nf-core/cutandrun documentation](https://nf-co.re/cutandrun/3.2.1) for further details regarding pipeline installation, operation, and parameter options.


# Conclusion 

Here we have presented nf-core/cutandrun, a best-practice bioinformatics pipeline for the end-to-end analysis of CUT&RUN, CUT&Tag and TIP-seq data. By adhering to nf-core community guidelines, the pipeline fulfils our defined prerequisites for modern bioinformatics analysis, and supports a range of execution platforms, operating systems, and resource-limited scenarios. Our pipeline supports user-defined spike-in genomes, offers additional read count-based normalisation modes, supports three distinct peak calling methods, and uniquely facilitates the analysis of TIP-seq data. We provide additional command line parameters for workflow customisation, including the specification of technical replicates, flow-switching options, and pipeline execution reporting. Collectively, we conclude that these properties of nf-core/cutandrun provide significant advantages relative to current alternative pipelines. Readers are encouraged to engage with the nf-core community through the official nf-core GitHub repository and Slack channels.


# Acknowledgements

This work was supported by the Francis Crick Institute, which receives its core funding from Cancer Research UK (CC001051, CC0102, FC010110), the UK Medical Research Council (CC001051, CC0102, FC010110), and the Wellcome Trust (CC001051, CC0102, FC010110). It was also supported by the Wellcome Trust (215593/Z/19/Z to J.U. and N.M.L.). N.M.L. was a Winton Group Leader in recognition of the Winton Charitable Foundation’s support towards the establishment of the Francis Crick Institute. For the purpose of Open Access, the author has applied a CC BY public copyright licence to any Author Accepted Manuscript version arising from this submission. We would like to thank all members of the Francis Crick Institute and the nf-core community for using the pipeline and providing structured feedback throughout its development and release.

# References