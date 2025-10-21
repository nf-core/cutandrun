#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2023 @chris-cheshire
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Author: @chris-cheshire

import os
import errno
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Create IGV session file from a list of files and associated colours - ".bed", ".bw", ".bigwig", ".tdf", ".gtf" files currently supported.'
Epilog = """Example usage: python igv_files_to_session.py <XML_OUT> <LIST_FILE> <GENOME>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument("XML_OUT", help="XML output file.")
argParser.add_argument(
    "LIST_FILE", help="Tab-delimited file containing two columns i.e. file_name\tcolour. Header isnt required."
)
argParser.add_argument(
    "GENOME", help="Full path to genome fasta file or shorthand for genome available in IGV e.g. hg19."
)
argParser.add_argument("GTF_BED")

## OPTIONAL PARAMETERS
argParser.add_argument(
    "-pp",
    "--path_prefix",
    type=str,
    dest="PATH_PREFIX",
    default="",
    help="Path prefix to be added at beginning of all files in input list file.",
)
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


def igv_files_to_session(XMLOut, ListFile, Genome, GtfBed, PathPrefix=""):
    makedir(os.path.dirname(XMLOut))

    fileList = []
    fin = open(ListFile, "r")
    while True:
        line = fin.readline()
        if line:
            ifile, colour = line.strip().split("\t")
            if len(colour.strip()) == 0:
                colour = "0,0,178"
            fileList.append((PathPrefix.strip() + ifile, colour))
        else:
            break
            fout.close()

    ## Construct groups
    groups = {}
    group_num = 1
    for ifile, colour in fileList:
        group = os.path.basename(ifile).split("_R")[0]
        if group not in groups:
            groups[group] = group_num
            group_num = group_num + 1

    ## ADD RESOURCES SECTION
    XMLStr = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    XMLStr += '<Session genome="%s" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">\n' % (Genome)
    XMLStr += "\t<Resources>\n"
    for ifile, colour in fileList:
        XMLStr += '\t\t<Resource path="%s"/>\n' % (ifile)
    XMLStr += '\t\t<Resource path="./%s"/>\n' % (GtfBed)
    XMLStr += "\t</Resources>\n"

    XMLStr += '\t<Panel height="1160" name="DataPanel" width="1897">\n'

    # Render gene file first
    XMLStr += (
        '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,48,73" '
    )
    XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" '
    XMLStr += (
        'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
        % (GtfBed, os.path.basename(GtfBed))
    )

    ## Do a GTF pass first
    for ifile, colour in fileList:
        extension = os.path.splitext(ifile)[1].lower()
        if extension in [".gtf"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )
        elif extension in [".gff"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )

    ## Then beds/narrowpeak
    for ifile, colour in fileList:
        extension = os.path.splitext(ifile)[1].lower()
        if extension in [".bed", ".broadpeak", ".narrowpeak"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="12" height="20" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )
        elif extension in [".bw", ".bigwig", ".tdf", ".bedGraph", ".bedgraph"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="true" autoscaleGroup="%s" clazz="org.broad.igv.track.DataSourceTrack" color="%s" '
                % (groups[os.path.basename(ifile).split("_R")[0]], colour)
            )
            XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="100" '
            XMLStr += (
                'id="%s" name="%s" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean">\n'
                % (ifile, os.path.basename(ifile))
            )
            XMLStr += '\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10" minimum="0.0" type="LINEAR"/>\n'
            XMLStr += "\t\t</Track>\n"
        elif extension in [".bam"]:
            pass
        # else:
        #     XMLStr += '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" ' % (colour)
        #     XMLStr += 'displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="10" height="20" '
        #     XMLStr += 'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n' % (ifile,os.path.basename(ifile))

    XMLStr += "\t</Panel>\n"
    # XMLStr += '\t<HiddenAttributes>\n\t\t<Attribute name="DATA FILE"/>\n\t\t<Attribute name="DATA TYPE"/>\n\t\t<Attribute name="NAME"/>\n\t</HiddenAttributes>\n'
    XMLStr += "</Session>"
    XMLOut = open(XMLOut, "w")
    XMLOut.write(XMLStr)
    XMLOut.close()


############################################
############################################
## RUN FUNCTION
############################################
############################################

igv_files_to_session(
    XMLOut=args.XML_OUT, ListFile=args.LIST_FILE, Genome=args.GENOME, GtfBed=args.GTF_BED, PathPrefix=args.PATH_PREFIX
)

############################################
############################################
############################################
############################################
