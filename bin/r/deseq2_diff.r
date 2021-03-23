#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## TODO

## TODO: TIDY CAMEL CASE AND SNAKE CASE

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(stringr)
library(BiocParallel)
library(magrittr)
library(GenomicRanges)
library(chromVAR)
library(DESeq2)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-c", "--control"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
    make_option(c("-t", "--treatment"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
    make_option(c("-b", "--bed"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
    make_option(c("-a", "--bam"    ), type="character", default=NULL    , metavar="path"   , help="TODO"  ),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."  ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."     ),
    make_option(c("-s", "--count_thresh"     ), type="integer", default='5', metavar="string" , help="TODO"     ),
    make_option(c("-@", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."   )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$control)){
    print_help(opt_parser)
    stop("Please provide a control group name.", call.=FALSE)
}

if (is.null(opt$treatment)){
    print_help(opt_parser)
    stop("Please provide a treatment group name.", call.=FALSE)
}

if (is.null(opt$bed)){
    print_help(opt_parser)
    stop("Please provide a list of bed files to load.", call.=FALSE)
}

if (is.null(opt$bam)){
    print_help(opt_parser)
    stop("Please provide a list of bam files to load", call.=FALSE)
}

# Get file lists
bed_list <- unlist(strsplit(opt$bed, ","))
bam_list <- unlist(strsplit(opt$bam, ","))

# Check same length
if (length(bed_list) != length(bam_list)) {
    print_help(opt_parser)
    stop("Bed and bam file list are different lengths", call.=FALSE)
}

################################################
################################################
## READ IN DATA FILES                         ##
################################################
################################################

# Create group list
groups = c(opt$control, opt$treatment)

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

################################################
################################################
## READ IN DATA FILES                         ##
################################################
################################################

# Init
mPeak = GRanges()
file_count = 0
file_list=vector()
# Read in bed files that match the control or treatment group
for(group in groups){
    search_res <-  str_detect(bed_list, group)
    file_list <- bed_list[search_res] %>% append(file_list, .)
}    
for(file in file_list) {
    peakRes = read.table(file, header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
    file_count = file_count + 1
}

# Create replicate counts and names
group_count = length(groups)
rep_count = file_count / group_count
reps = paste0("rep", 1:rep_count)

# Create peak table and count matrix
masterPeak = reduce(mPeak)
countMat = matrix(NA, length(masterPeak), file_count)
colnames(countMat) = paste(rep(groups, 2), rep(reps, each = 2), sep = "_")

# Read in bam files that match the control or treatment group
for(i in seq_along(groups)){
    search_res <-  str_detect(bam_list, groups[i])
    file_list <- bam_list[search_res]
    print(file_list)
    
    for(j in seq_along(file_list)) {
        fragment_counts <- getCounts(file_list[j], masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
        countMat[, (((i-1)*group_count) + (j-1)) + 1] = counts(fragment_counts)[,1]
    }
}


################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

selectR = which(rowSums(countMat) > opt$count_thresh) ## Create index list for peak count filter
dataS = countMat[selectR,] ## Select data from filter
condition = factor(rep(groups, each = length(reps)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


################################################
################################################
## FORMAT RESULTS                             ##
################################################
################################################

countMatDiff = cbind(dataS, normDDS, res) ## Combine deseq results
countMatDiff

peakFiltered = masterPeak[selectR,]
peakFiltered

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

if(FALSE) {
################################################
################################################
## OUTPUT RESULTS                             ##
################################################
################################################

#write.csv(countMatDiff,"count_mat_diff.csv", row.names=FALSE)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}
}
################################################
################################################
################################################
################################################