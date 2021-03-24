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
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

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
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog." ),
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

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}
setwd(opt$outdir)

## Create index list for peak count filter
selectR = which(rowSums(countMat) > opt$count_thresh) 
dataS = countMat[selectR,] ## Select data from filter
condition = factor(rep(groups, each = length(reps)))

samples.vec <- sort(colnames(countMat))
groups      <- sub("_[^_]+$", "", samples.vec)
if (length(unique(groups)) == 1 || length(unique(groups)) == length(samples.vec)) {
    quit(save = "no", status = 0, runLast = FALSE)
}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")
if (file.exists(DDSFile) == FALSE) {
    dds = DESeqDataSetFromMatrix(countData = dataS,colData = DataFrame(condition),design = ~ condition)
    dds     <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
    if (!opt$vst) {
        vst_name <- "rlog"
        rld      <- rlog(dds)
    } else {
        vst_name <- "vst"
        rld      <- varianceStabilizingTransformation(dds)
    }
    assay(dds, vst_name) <- assay(rld)
    save(dds,file=DDSFile)
} else {
    load(DDSFile)
    vst_name <- intersect(assayNames(dds), c("vst", "rlog"))
    if (length(vst_name)==0) { # legacy might mean vst was saved as a separate object called rld
        vst_name <- "loaded_rld"
        assay(dds, vst_name) <- assay(rld)
    } else {
        vst_name==vst_name[1]
    }
}

## Normalised results
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
countMatDiff = cbind(dataS, normDDS, res) ## Combine deseq results
write.csv(countMatDiff, paste(opt$outprefix, "count_mat_diff.csv", sep=""), row.names=FALSE)

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on, and adds an assessment of how well
##' each PC explains the experimental grouping of the data.
##' 
##' @param object The DESeq2DataSet object.
##' @param intgroup interesting groups: a character vector of names in 'colData(x)' to use for grouping.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object, intgroup = "condition", ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }  else {
        colData(object)[[intgroup]]
    }
    d <- cbind(pca$x, group = group, intgroup.df, name = colnames(object))
    percentFrame <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar, groupR=0.0)
    for (ipc in seq(along=percentVar)) {
        fit1 <- lm(pca$x[,ipc]  ~ group)
        percentFrame$groupR[ipc] <- 100*summary(fit1)$r.squared
    }
    attr(d, "percentVar") <- percentFrame
    return(d)
}

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")
if (file.exists(PlotFile) == FALSE) {
    pdf(file=PlotFile,onefile=TRUE,width=7,height=7)
    
    ## PCA
    ntop <- c(500, Inf)
    for (n_top_var in ntop) {
        pca.data      <- plotPCA_vst(dds, assay=vst_name,intgroup=c("condition"),ntop=n_top_var)
        percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
        plot_subtitle <- ifelse(n_top_var==Inf, "All peaks", paste("Top", n_top_var, "peaks"))
        pl <- ggplot(pca.data, aes(PC1, PC2, color=condition)) +
            geom_point(size=3) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) +
            labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) + 
            theme(legend.position="top",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
        print(pl)
        
        pl <- ggplot(attr(pca.data, "percentVar"), aes(x=PC, y=percentVar)) +
            geom_line(aes(colour="explained by PC")) +
            geom_line(aes(y=groupR, colour="of PC explained by condition")) +
            scale_x_continuous(breaks=seq(along=percentVar), minor_breaks=NULL)  +
            labs(title="Diagnostics of PCs", subtitle=plot_subtitle, x="Component", y="Percentage explaned", colour="Percentage variation") +
            theme_bw() +
            theme(legend.position="top")
        print(pl)
        
        pc_r <- order(attr(pca.data, "percentVar")$groupR, decreasing=TRUE)
        pl <- ggplot(pca.data, aes_string(paste0("PC", pc_r[1]), paste0("PC", pc_r[2]), color="condition")) +
            geom_point(size=3) +
            xlab(paste0("PC", pc_r[1], ": ",percentVar[pc_r[1]],"% variance")) +
            ylab(paste0("PC", pc_r[2], ": ",percentVar[pc_r[2]],"% variance")) +
            labs(title = paste0("Group-Explanatory PCs of ", vst_name, "-tranformed data"), subtitle = plot_subtitle) + 
            theme(legend.position="top",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
        print(pl)
    } # at end of loop, we'll be using the user-defined ntop if any, else all peaks
    
    ## VOLCANO PLOT
    plotMA(dds)
    
    ## WRITE PC1 vs PC2 VALUES TO FILE
    pca.vals           <- pca.data[,1:2]
    colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
    pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
    write.table(pca.vals,file=paste(opt$outprefix,".pca.vals.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)
    
    ## SAMPLE CORRELATION HEATMAP
    sampleDists      <- dist(t(assay(dds, vst_name)))
    sampleDistMatrix <- as.matrix(sampleDists)
    colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(
        sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colors,
        main=paste("Euclidean distance between", vst_name, "of samples")
    )
    
    ## WRITE SAMPLE DISTANCES TO FILE
    write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix,".sample.dists.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    
    dev.off()
}

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################
#if(FALSE){
SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir,recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir,opt$outprefix,".size_factors.RData",sep="")
if (file.exists(NormFactorsFile) == FALSE) {
    normFactors <- sizeFactors(dds)
    save(normFactors,file=NormFactorsFile)
    
    for (name in names(sizeFactors(dds))) {
        sizeFactorFile <- paste(SizeFactorsDir,name,".txt",sep="")
        if (file.exists(sizeFactorFile) == FALSE) {
            write(as.numeric(sizeFactors(dds)[name]),file=sizeFactorFile)
        }
    }
}
#}

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

################################################
################################################
################################################
################################################