#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/design.csv


import os
import sys
from collections import Counter
from pathlib import Path


logger = logging.getLogger()

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("IGG", help="Boolean for whether or not igg is given")
    return parser.parse_args(args)

class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "At least the first FASTQ file is required."
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            assert (
                Path(row[self._first_col]).suffixes == Path(row[self._second_col]).suffixes
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

def check_samplesheet(file_in, file_out, igg_control):
    """
    Detect the tabular format.

    group,replicate,control_group,fastq_1,fastq_2
    WT,1,1,WT_LIB1_REP1_1.fastq.gz,WT_LIB1_REP1_2.fastq.gz
    WT,1,1,WT_LIB2_REP1_1.fastq.gz,WT_LIB2_REP1_2.fastq.gz
    WT,2,1,WT_LIB1_REP2_1.fastq.gz,WT_LIB1_REP2_2.fastq.gz
    KO,1,2,KO_LIB1_REP1_1.fastq.gz,KO_LIB1_REP1_2.fastq.gz
    IGG,1,1,KO_LIB1_REP1_1.fastq.gz,IGG_LIB1_REP1_2.fastq.gz
    IGG,2,2,KO_LIB1_REP1_1.fastq.gz,IGG_LIB1_REP1_2.fastq.gz
    """
    peek = handle.read(2048)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    handle.seek(0)
    return dialect

    igg_present = False

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = ["group", "replicate", "control_group", "fastq_1", "fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

            if not igg_present:
                if 'igg' in lspl:
                    igg_present = True

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, replicate, control_group, fastq_1, fastq_2 = lspl[: len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            ## Check replicate entry is integer
            if not replicate.isdigit():
                print_error("Replicate id not an integer!", "Line", line)
            replicate = int(replicate)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = [sample, str(replicate), control_group, "0", fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = [sample, str(replicate), control_group, "1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)
            ## Create sample mapping dictionary = {sample: {replicate : [ single_end, fastq_1, fastq_2 ]}}
            if sample not in sample_run_dict:
                sample_run_dict[sample] = {}
            if replicate not in sample_run_dict[sample]:
                sample_run_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_run_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_run_dict[sample][replicate].append(sample_info)

    ## Check igg_control parameter is consistent with input groups
    if (igg_control == 'true' and not igg_present):
        print("ERROR: No 'igg' group was found in " + str(file_in) + " If you are not supplying an IgG control, please specify --igg_control 'false' on command line.")
        sys.exit(1)

    if (igg_control == 'false' and igg_present):
        print("ERROR: Parameter --igg_control was set to false, but an 'igg' group was found in " + str(file_in) + ".")
        sys.exit(1)

    ## Check control groups have unique ids that are the same as their replicate ids
    control_group_ids = []
    if igg_present:
        for key, data in sample_run_dict["igg"].items():
            for tech_rep in data:
                if tech_rep[2] not in control_group_ids:
                    control_group_ids.append(tech_rep[2])

                if(tech_rep[2] != str(key)):
                    print("ERROR: IgG groups must have a control id equal to the replicate id")
                    sys.exit(1)

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:

            fout.write(",".join(["id", "group", "replicate", "control_group", "single_end", "fastq_1", "fastq_2"]) + "\n")
            for sample in sorted(sample_run_dict.keys()):

                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_run_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1..<num_replicates>!",
                        "Group",
                        sample,
                    )
                for replicate in sorted(sample_run_dict[sample].keys()):

                    ## Check control group exists
                    if igg_present:
                        for tech_rep in sample_run_dict[sample][replicate]:
                            if tech_rep[2] not in control_group_ids:
                                print_error(
                                    "Control group does not exist",
                                    tech_rep[2]
                                )

                        ## Check tech reps have same control group id
                        check_group = sample_run_dict[sample][replicate][0][2]
                        for tech_rep in sample_run_dict[sample][replicate]:
                            if tech_rep[2] != check_group:
                                print_error(
                                    "Control group must match within technical replicates",
                                    tech_rep[2]
                                )

                    ## Check that multiple runs of the same sample are of the same datatype
                    if not all(
                        x[0] == sample_run_dict[sample][replicate][0][0] for x in sample_run_dict[sample][replicate]
                    ):
                        print_error(
                            "Multiple runs of a sample must be of the same datatype!",
                            "Group",
                            sample,
                        )
                    ## Write to file
                    for idx, sample_info in enumerate(sample_run_dict[sample][replicate]):
                        sample_id = "{}_R{}_T{}".format(sample, replicate, idx + 1)
                        fout.write(",".join([sample_id] + sample_info) + "\n")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT, args.IGG)

if __name__ == "__main__":
    sys.exit(main())
