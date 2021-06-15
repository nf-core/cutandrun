#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## DIFFERENTIAL ANALYSIS, SCATTERPLOTS, HEATMAP AND PCA FROM COUNT FILE
    ## - COUNT FILE IS PRODUCED FROM SEACR PEAK BED FILES AND ALIGNMENT BAMS
    ## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R

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
library(lattice)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-g", "--groups"), type="character", default=NULL    , metavar="string"   , help="comma-separated list of experimental group names" ),
    make_option(c("-b", "--bed"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
    make_option(c("-a", "--bam"    ), type="character", default=NULL    , metavar="path"   , help="TODO"  ),
    make_option(c("-i", "--include"), type="character", default=NULL    , metavar="string"   , help="experimental groups to include in analysis"),
    make_option(c("-e", "--exclude"), type="character", default=NULL    , metavar="string"   , help="experimental groups to exclude in analysis"),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."  ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."     ),
    make_option(c("-s", "--count_thresh"     ), type="integer", default='1', metavar="string" , help="TODO"     ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog." ),
    make_option(c("-@", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."   )
)

#make_option(c("-c", "--control"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
#make_option(c("-t", "--treatment"    ), type="character", default=NULL    , metavar="path"   , help="TODO" ),
opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$groups)){
    print_help(opt_parser)
    stop("Please provide group names for analysis.", call.=FALSE)
}

#if (is.null(opt$control)){
#    print_help(opt_parser)
#    stop("Please provide a control group name.", call.=FALSE)
#}
#
#if (is.null(opt$treatment)){
#    print_help(opt_parser)
#    stop("Please provide a treatment group name.", call.=FALSE)
#}

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

# Check same length <- MAYBE TAKE THIS OUT IF WE SIMPLY PASS THE ALIGNMENT CHANNEL AS IT WILL CONTAIN IGG
#if (length(bed_list) != length(bam_list)) {
#    print_help(opt_parser)
#    stop("Bed and bam file list are different lengths", call.=FALSE)
#}

################################################
################################################
## READ IN DATA FILES                         ##
################################################
################################################

# Create group list
# groups = c(opt$control, opt$treatment) <- groups is now just opt$group
groups = unlist(strsplit(opt$groups, split=","))
if (!is.null(opt$include)) {
    groups = unlist(strsplit(opt$include, split=","))
}
if (!is.null(opt$exclude)) {
    exclude_vec = unlist(strsplit(opt$exclude, split=","))
    matching = match(groups, exclude_vec)
    groups = groups[!matching]
}

################ redo from here #######################
#print(groups)
#print(bed_list)

# Init
mPeak = GRanges()
file_count = 0
file_list=vector()
sample_mat = matrix(NA, length(bam_list), 2)
# Read in bed files that match the control or treatment group
for(group in groups){
    search_res <-  str_detect(bed_list, group)
    if (!any(search_res)) {
        stop(paste("group", group, "was not found amongst bed files"))
    }
    file_list <- bed_list[search_res] %>% append(file_list, .)
}
k=1
for(file in file_list) {
    peakRes = read.table(file, header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
    file_count = file_count + 1
    file_name = basename(file)
    file_name = strsplit(file_name, "[.]")[[1]][1]
    file_split = strsplit(file_name, "[_]")[[1]]
    file_group = paste(file_split[-length(file_split)], collapse = "_")
    file_rep = file_split[length(file_split)]
    file_rep_num = gsub("[^0-9.]", "",  file_rep)
    sample_mat[k,1] = file_group
    sample_mat[k,2] = file_rep_num
    k=k+1
}


# Create peak table and count matrix
masterPeak = reduce(mPeak)
countMat = matrix(NA, length(masterPeak), file_count)
sample_list = paste(sample_mat[,1], sample_mat[,2], sep = "_R")
colnames(countMat) = sample_list
# Read in bam files that match the control or treatment group
for(i in seq_along(groups)){
    search_res <-  str_detect(bam_list, groups[i])
    if (!any(search_res)) {
        stop(paste("group", i, "was not found amongst bam files"))
    }
}

for(j in seq_along(bam_list)) {
    rep_search <- str_detect(bam_list, sample_list[j])
    file_now <- bam_list[rep_search]
    fragment_counts <- getCounts(file_now, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[,j] = counts(fragment_counts)[,1]
}

#################################################################

# # Init
# mPeak = GRanges()
# file_count = 0
# file_list=vector()
# # Read in bed files that match the control or treatment group
# for(group in groups){
#     search_res <-  str_detect(bed_list, group)
#     if (!any(search_res)) {
#         stop(paste("group", group, "was not found amongst bed files"))
#     }
#     file_list <- bed_list[search_res] %>% append(file_list, .)
# }
# for(file in file_list) {
#     peakRes = read.table(file, header = FALSE, fill = TRUE)
#     mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
#     file_count = file_count + 1
# }

###########################################################################




# # Create replicate counts and names
# group_count = length(groups)
# rep_count = file_count / group_count
# reps = paste0("rep", 1:rep_count)
# expected_rep_str = paste0("_R", 1:rep_count)

# # Create peak table and count matrix
# masterPeak = reduce(mPeak)
# countMat = matrix(NA, length(masterPeak), file_count)
# colnames(countMat) = paste(rep(groups, rep_count), rep(reps, each = group_count, sep = "_"))

# # Read in bam files that match the control or treatment group
# for(i in seq_along(groups)){
#     search_res <-  str_detect(bam_list, groups[i])
#     if (!any(search_res)) {
#         stop(paste("group", i, "was not found amongst bam files"))
#     }
#     file_list <- bam_list[search_res]

#     for(j in seq_along(file_list)) {
#         rep_search <- str_detect(file_list, expected_rep_str[j])
#         file_now <- file_list[rep_search]
#         fragment_counts <- getCounts(file_now, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
#         countMat[, ((j*group_count) - (group_count - i))] = counts(fragment_counts)[,1]
#     }
# }

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
condition = factor(sample_mat[,1])

samples.vec <- sort(colnames(countMat))

#groups      <- sub("_[^_]+$", "", samples.vec)
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
## NEW PLOT QC                                ##
################################################
################################################
# if (FALSE) {
# PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")
# if (file.exists(PlotFile) == FALSE) {
#     pdf(file=PlotFile,onefile=TRUE,width=7,height=7)

#     ## PCA
#     pca.data <- DESeq2::plotPCA(rld,intgroup=c("condition"),returnData=TRUE)
#     percentVar <- round(100 * attr(pca.data, "percentVar"))
#     plot <- ggplot(pca.data, aes(PC1, PC2, color=condition)) +
#         geom_point(size=3) +
#         xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#         ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#         theme(panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 panel.background = element_blank(),
#                 panel.border = element_rect(colour = "black", fill=NA, size=1))
#     print(plot)

#     ## WRITE PC1 vs PC2 VALUES TO FILE
#     pca.vals <- pca.data[,1:2]
#     colnames(pca.vals) <- paste(colnames(pca.vals),paste(percentVar,'% variance',sep=""), sep=": ")
#     pca.vals <- cbind(sample = rownames(pca.vals), pca.vals)
#     write.table(pca.vals,file=paste(opt$outprefix,".pca.vals.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=TRUE)

#     ## SAMPLE CORRELATION HEATMAP
#     sampleDists <- dist(t(assay(rld)))
#     sampleDistMatrix <- as.matrix(sampleDists)
#     colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#     pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)

#     ## WRITE SAMPLE DISTANCES TO FILE
#     write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix,".sample.dists.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

#     dev.off()
# }
# }

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
        pca.data      <- plotPCA_vst(dds, assay=vst_name,intgroup="condition",ntop=n_top_var)
        percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
        plot_subtitle <- ifelse(n_top_var==Inf, "All peaks", paste("Top", n_top_var, "peaks"))
        # PCA PLOT 1 - FIRST PC
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

        # PCA PLOT 2 - DIAGNOSTIC OF PCs
        pl <- ggplot(attr(pca.data, "percentVar"), aes(x=PC, y=percentVar)) +
            geom_line(aes(colour="explained by PC")) +
            geom_line(aes(y=groupR, colour="of PC explained by condition")) +
            scale_x_continuous(breaks=seq(along=percentVar), minor_breaks=NULL)  +
            labs(title="Diagnostics of PCs", subtitle=plot_subtitle, x="Component", y="Percentage explaned", colour="Percentage variation") +
            theme_bw() +
            theme(legend.position="top")
        print(pl)

        # PCA PLOT 3 - GROUP-EXPLANATORY PCs
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
## LOOP THROUGH COMPARISONS                   ##
################################################
################################################

# if (FALSE) {
#     ResultsFile <- paste(opt$outprefix,".results.txt",sep="")
#     if (file.exists(ResultsFile) == FALSE) {
#         raw.counts <- counts(dds,normalized=FALSE)
#         colnames(raw.counts) <- paste(colnames(raw.counts),'raw',sep='.')
#         pseudo.counts <- counts(dds,normalized=TRUE)
#         colnames(pseudo.counts) <- paste(colnames(pseudo.counts),'pseudo',sep='.')

#         deseq2_results_list <- list()
#         comparisons <- combn(unique(groups),2)
#         for (idx in 1:ncol(comparisons)) {
#             control.group <- comparisons[1,idx]
#             treat.group <- comparisons[2,idx]
#             CompPrefix <- paste(control.group,treat.group,sep="vs")
#             cat("Saving results for ",CompPrefix," ...\n",sep="")

#             CompOutDir <- paste(CompPrefix,'/',sep="")
#             if (file.exists(CompOutDir) == FALSE) {
#                 dir.create(CompOutDir,recursive=TRUE)
#             }

#             control.samples <- samples.vec[which(groups == control.group)]
#             treat.samples <- samples.vec[which(groups == treat.group)]
#             comp.samples <- c(control.samples,treat.samples)

#             comp.results <- results(dds,contrast=c("condition",c(control.group,treat.group)))
#             comp.df <- as.data.frame(comp.results)
#             #comp.table <- cbind(interval.table, as.data.frame(comp.df), raw.counts[,paste(comp.samples,'raw',sep='.')], pseudo.counts[,paste(comp.samples,'pseudo',sep='.')])
#             comp.table <- cbind(as.data.frame(comp.df), raw.counts[,paste(comp.samples,'raw',sep='.')], pseudo.counts[,paste(comp.samples,'pseudo',sep='.')])

#             ## WRITE RESULTS FILE
#             CompResultsFile <- paste(CompOutDir,CompPrefix,opt$outsuffix,".deseq2.results.txt",sep="")
#             write.table(comp.table, file=CompResultsFile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

#             ## FILTER RESULTS BY FDR & LOGFC AND WRITE RESULTS FILE
#             pdf(file=paste(CompOutDir,CompPrefix,opt$outsuffix,".deseq2.plots.pdf",sep=""),width=10,height=8)
#             if (length(comp.samples) > 2) {
#                 for (MIN_FDR in c(0.01,0.05)) {

#                     ## SUBSET RESULTS BY FDR
#                     #pass.fdr.table <- subset(comp.table, padj < MIN_FDR)
#                     #pass.fdr.up.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange > 0)
#                     #pass.fdr.down.table <- subset(comp.table, padj < MIN_FDR & log2FoldChange < 0)

#                     ## SUBSET RESULTS BY FDR AND LOGFC
#                     #pass.fdr.logFC.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1)
#                     #pass.fdr.logFC.up.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange > 0)
#                     #pass.fdr.logFC.down.table <- subset(comp.table, padj < MIN_FDR & abs(log2FoldChange) >= 1 & log2FoldChange < 0)

#                     ## WRITE RESULTS FILE
#                     #CompResultsFile <- paste(CompOutDir,CompPrefix,opt$outsuffix,".deseq2.FDR",MIN_FDR,".results.txt",sep="")
#                     #CompBEDFile <- paste(CompOutDir,CompPrefix,opt$outsuffix,".deseq2.FDR",MIN_FDR,".results.bed",sep="")
#                     #write.table(pass.fdr.table, file=CompResultsFile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
#                     #write.table(pass.fdr.table[,c("Chr","Start","End","Geneid","log2FoldChange","Strand")], file=CompBEDFile, col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)

#                     ## MA PLOT & VOLCANO PLOT
#                     DESeq2::plotMA(comp.results, main=paste("MA plot FDR <= ",MIN_FDR,sep=""), ylim=c(-2,2),alpha=MIN_FDR)
#                     #plot(comp.table$log2FoldChange, -1*log10(comp.table$padj), col=ifelse(comp.table$padj<=MIN_FDR, "red", "black"), xlab="logFC", ylab="-1*log10(FDR)", main=paste("Volcano plot FDR <=",MIN_FDR,sep=" "), pch=20)

#                     ## ADD COUNTS TO LOGFILE
#                     #cat(CompPrefix," genes with FDR <= ",MIN_FDR,": ",nrow(pass.fdr.table)," (up=",nrow(pass.fdr.up.table),", down=",nrow(pass.fdr.down.table),")","\n",file=LogFile,append=TRUE,sep="")
#                     #cat(CompPrefix," genes with FDR <= ",MIN_FDR," & FC > 2: ",nrow(pass.fdr.logFC.table)," (up=",nrow(pass.fdr.logFC.up.table),", down=",nrow(pass.fdr.logFC.down.table),")","\n",file=LogFile,append=TRUE,sep="")

#                 }
#                 cat("\n",file=LogFile,append=TRUE,sep="")
#             }

#             ## SAMPLE CORRELATION HEATMAP
#             rld.subset <- assay(rld)[,comp.samples]
#             sampleDists <- dist(t(rld.subset))
#             sampleDistMatrix <- as.matrix(sampleDists)
#             colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#             pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)

#             ## SCATTER PLOT FOR RLOG COUNTS
#             combs <- combn(comp.samples,2,simplify=FALSE)
#             clabels <- sapply(combs,function(x){paste(x,collapse=' & ')})
#             plotdat <- data.frame(x=unlist(lapply(combs, function(x){rld.subset[, x[1] ]})),y=unlist(lapply(combs, function(y){rld.subset[, y[2] ]})),comp=rep(clabels, each=nrow(rld.subset)))
#             plot <- xyplot(y~x|comp,plotdat,
#                     panel=function(...){
#                         panel.xyplot(...)
#                         panel.abline(0,1,col="red")
#                     },
#                     par.strip.text=list(cex=0.5))
#             print(plot)
#             dev.off()

#             colnames(comp.df) <- paste(CompPrefix,".",colnames(comp.df),sep="")
#             deseq2_results_list[[idx]] <- comp.df
#         }

#         ## WRITE RESULTS FROM ALL COMPARISONS TO FILE
#         deseq2_results_table <- cbind(do.call(cbind, deseq2_results_list),raw.counts,pseudo.counts)
#         write.table(deseq2_results_table, file=ResultsFile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
#     }
# }

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

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
#}
################################################
################################################
################################################
################################################
