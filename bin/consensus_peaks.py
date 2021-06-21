#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on June 29th 2018 to annotate merged peaks
#######################################################################
#######################################################################

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import argparse
import upsetplot

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = 'Upset ven diagram of consensus peaks.'
Epilog = """Example usage: python consensus_peaks.py <MERGED_INTERVAL_FILE> """

parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
parser.add_argument('--peaks', help="Merged peaks interval file with replicate counts column.")
# parser.add_argument('--output', help="Full path to output directory.")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# create list of data frames, one for each group consensus peaks file
peak_file_list = glob.glob(args.peaks)
peak_df_list = list()
for i in list(range(len(peak_file_list))):
    peaks_i = pd.read_csv(peak_file_list[i], sep='\t', header=None, usecols=[0,1,2,8,9], names=['chrom','start','end','sample_reps','count'])
    peaks_i['sample_reps'] = peaks_i['sample_reps'].replace(".peaks.bed.stringent.bed", "", regex=True)
    peak_df_list.append(peaks_i)
    # print("peak i")
    # print(peaks_i.head(10))

    reps2 = peaks_i[peaks_i["count"]>1]
    # print("reps2")
    # print(reps2.head(10))

# print(peak_df_list)

# add sorted column to each dataframe, and make new condensed dataframe list
summary_peak_df_list = list()
for i in list(range(len(peak_df_list))):
    peaks_i = peak_df_list[i]
    peaks_i['sorted_samples'] = ''
    rows_now = peaks_i.shape[0]
    for j in list(range(rows_now)):
        sample_list = peaks_i.at[j,'sample_reps']
        sample_array = np.unique(sample_list.split(','))
        sample_sorted = sorted(sample_array)
        sample_str = ",".join(sample_sorted)
        peaks_i.at[j,'sorted_samples'] = sample_str
    summary_peaks_i = peaks_i[['sorted_samples', 'count']].groupby(['sorted_samples'], as_index = False).sum()
    summary_peak_df_list.append(summary_peaks_i)
    # print(summary_peaks_i)
# print(summary_peak_df_list)


# Construct data in appropriate format for upsetplot, and plot
for i in list(range(len(summary_peak_df_list))):
    df_i = summary_peak_df_list[i]
    print(type(df_i))
    print(df_i)
    categories = df_i.shape[0]
    cat_list = np.zeros(categories)
    for j in list(range(categories)):
        cat_list[j] = df_i.at[j,'sorted_samples'].split(',')
        # print(blop)
    print(cat_list)

    # Plot
    # peak_counts = upsetplot.from_memberships(cat_list, data = df_i['count'])
    # upsetplot.plot(peak_counts)
    # plt.show()



if (False):
    ############################################
    ############################################
    ## PARSE ARGUMENTS
    ############################################
    ############################################
    Description = 'Aggregate columns from merged peak file, filter for minimum replicate threshold, and produce Upset ven diagram of consensus peaks.'
    Epilog = """Example usage: python threshold_peaks.py <MERGED_INTERVAL_FILE> <SAMPLE_NAME_LIST> <OUTFILE> --min_replicates 1"""

    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    ## REQUIRED PARAMETERS
    argParser.add_argument('MERGED_INTERVAL_FILE', help="Merged MACS2 interval file created using linux sort and mergeBed.")
    argParser.add_argument('SAMPLE_NAME_LIST', help="Comma-separated list of sample names as named in individual MACS2 broadPeak/narrowPeak output file e.g. SAMPLE_R1 for SAMPLE_R1_peak_1.")
    argParser.add_argument('OUTFILE', help="Full path to output directory.")

    ## OPTIONAL PARAMETERS
    argParser.add_argument('-in', '--is_narrow_peak', dest="IS_NARROW_PEAK", help="Whether merged interval file was generated from narrow or broad peak files (default: False).",action='store_true')
    argParser.add_argument('-mr', '--min_replicates', type=int, dest="MIN_REPLICATES", default=1, help="Minumum number of replicates per sample required to contribute to merged peak (default: 1).")
    args = argParser.parse_args()

    ############################################
    ############################################
    ## HELPER FUNCTIONS
    ############################################
    ############################################

    def makedir(path):

        if not len(path) == 0:
            try:
                os.makedirs(path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

    ############################################
    ############################################
    ## MAIN FUNCTION
    ############################################
    ############################################

    ## MergedIntervalTxtFile is file created using commands below:
    ## 1) broadPeak
    ## sort -k1,1 -k2,2n <MACS_BROADPEAK_FILES_LIST> | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt
    ## 2) narrowPeak
    ## sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt

    def macs2_merged_expand(MergedIntervalTxtFile,SampleNameList,OutFile,isNarrow=False,minReplicates=1):

        makedir(os.path.dirname(OutFile))

        combFreqDict = {}
        totalOutIntervals = 0
        SampleNameList = sorted(SampleNameList)
        fin = open(MergedIntervalTxtFile,'r')
        fout = open(OutFile,'w')
        oFields = ['chr','start','end','interval_id','num_peaks','num_samples'] + [x+'.bool' for x in SampleNameList] + [x+'.fc' for x in SampleNameList] + [x+'.qval' for x in SampleNameList] + [x+'.pval' for x in SampleNameList] + [x+'.start' for x in SampleNameList] + [x+'.end' for x in SampleNameList]
        if isNarrow:
            oFields += [x+'.summit' for x in SampleNameList]
        fout.write('\t'.join(oFields) + '\n')
        while True:
            line = fin.readline()
            if line:
                lspl = line.strip().split('\t')

                chromID = lspl[0]; mstart = int(lspl[1]); mend = int(lspl[2]);
                starts = [int(x) for x in lspl[3].split(',')]; ends = [int(x) for x in lspl[4].split(',')]
                names = lspl[5].split(','); fcs = [float(x) for x in lspl[8].split(',')]
                pvals = [float(x) for x in lspl[9].split(',')]; qvals = [float(x) for x in lspl[10].split(',')]
                summits = []
                if isNarrow:
                    summits = [int(x) for x in lspl[11].split(',')]

                ## GROUP SAMPLES BY REMOVING TRAILING *_R*
                groupDict = {}
                for sID in ['_'.join(x.split('_')[:-2]) for x in names]:
                    gID = '_'.join(sID.split('_')[:-1])
                    if gID not in groupDict:
                        groupDict[gID] = []
                    if sID not in groupDict[gID]:
                        groupDict[gID].append(sID)

                ## GET SAMPLES THAT PASS REPLICATE THRESHOLD
                passRepThreshList = []
                for gID,sIDs in groupDict.items():
                    if len(sIDs) >= minReplicates:
                        passRepThreshList += sIDs

                ## GET VALUES FROM INDIVIDUAL PEAK SETS
                fcDict = {}; qvalDict = {}; pvalDict = {}; startDict = {}; endDict = {}; summitDict = {}
                for idx in range(len(names)):
                    sample = '_'.join(names[idx].split('_')[:-2])
                    if sample in passRepThreshList:
                        if sample not in fcDict:
                            fcDict[sample] = []
                        fcDict[sample].append(str(fcs[idx]))
                        if sample not in qvalDict:
                            qvalDict[sample] = []
                        qvalDict[sample].append(str(qvals[idx]))
                        if sample not in pvalDict:
                            pvalDict[sample] = []
                        pvalDict[sample].append(str(pvals[idx]))
                        if sample not in startDict:
                            startDict[sample] = []
                        startDict[sample].append(str(starts[idx]))
                        if sample not in endDict:
                            endDict[sample] = []
                        endDict[sample].append(str(ends[idx]))
                        if isNarrow:
                            if sample not in summitDict:
                                summitDict[sample] = []
                            summitDict[sample].append(str(summits[idx]))

                samples = sorted(fcDict.keys())
                if samples != []:
                    numSamples = len(samples)
                    boolList  = ['TRUE' if x in samples else 'FALSE' for x in SampleNameList]
                    fcList = [';'.join(fcDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    qvalList = [';'.join(qvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    pvalList = [';'.join(pvalDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    startList = [';'.join(startDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    endList = [';'.join(endDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    oList = [str(x) for x in [chromID,mstart,mend,'Interval_'+str(totalOutIntervals+1),len(names),numSamples]+boolList+fcList+qvalList+pvalList+startList+endList]
                    if isNarrow:
                        oList += [';'.join(summitDict[x]) if x in samples else 'NA' for x in SampleNameList]
                    fout.write('\t'.join(oList) + '\n')

                    tsamples = tuple(sorted(samples))
                    if tsamples not in combFreqDict:
                        combFreqDict[tsamples] = 0
                    combFreqDict[tsamples] += 1
                    totalOutIntervals += 1

            else:
                fin.close()
                fout.close()
                break

        ## WRITE FILE FOR INTERVAL INTERSECT ACROSS SAMPLES.
        ## COMPATIBLE WITH UPSETR PACKAGE.
        fout = open(OutFile[:-4]+'.intersect.txt','w')
        combFreqItems = sorted([(combFreqDict[x],x) for x in combFreqDict.keys()],reverse=True)
        for k,v in combFreqItems:
            fout.write('%s\t%s\n' % ('&'.join(v),k))
        fout.close()

    ############################################
    ############################################
    ## RUN FUNCTION
    ############################################
    ############################################

    macs2_merged_expand(MergedIntervalTxtFile=args.MERGED_INTERVAL_FILE,SampleNameList=args.SAMPLE_NAME_LIST.split(','),OutFile=args.OUTFILE,isNarrow=args.IS_NARROW_PEAK,minReplicates=args.MIN_REPLICATES)

    ############################################
    ############################################
    ############################################
    ############################################
