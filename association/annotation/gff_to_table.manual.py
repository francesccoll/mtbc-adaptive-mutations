#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio import SeqFeature as sf
# from BCBio import GFF

# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------

# At the moment only GFF format is supported.
# Remove first "source" GFF line, e.g.: "MW2	EMBL/GenBank/SwissProt	source	1	2820462	"
# NOTE: gff_to_table.manual.py script has the same functionality as gff_to_table.py
# but does not make use of BCBio to parse the input GFF. Instead, gff_to_table.manual.py had to be
# written to deal with the format of Mycobacterium_tuberculosis_H37Rv_gff_v4.gff


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "Python script to parse GFF file to table format"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-i", "--input_gff", action="store", dest="input_gff",
        help="input genome annotation in GFF format",
        required=True, metavar="INPUT_GFF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output genome annotation in table format",
        required=True, metavar="OUTPUT_TABLE"
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------


def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Making sure input files exist
    if not os.path.isfile(args.input_gff):
        logging.error(f'Input genome {args.input_gff} not found!')
        sys.exit(-1)

    # Defining feature types and attributes to extract
    feature_types = ['CDS']
    # attributes = ['locus_tag', 'gene', 'product']
    attributes = ['Locus', 'Name', 'Product']

    # Creating table header
    header = attributes + ['type', 'start', 'end', 'strand', 'snpEff_locus_id']

    # Loading input genome
    logging.info(f"Opening input annotated genome {args.input_gff}")
    in_handle = open(args.input_gff)
    gff_tab_file = open(args.output_table, 'w')
    gff_tab_file.write('\t'.join(header)+'\n')
    # for rec in GFF.parse(in_handle):
    #     print(rec)
    #     for feature in rec.features:
    #         print(feature)
    for line in in_handle:
        print(line)
        items = line.strip().split('\t')
        feature_type = items[2]
        if feature_type in feature_types:
            chr_id = items[0]
            strand = items[6]
            locus_start = items[3]
            locus_end = items[4]
            ann = items[8]
            print('feature_type ' + feature_type + ' strand ' + strand + ' locus_start ' + locus_start +
                  ' locus_end ' + locus_end)
            ann_dict = dict()
            ann_items = ann.split(';')
            for i in ann_items:
                # NOTE, making symbol is present in item
                if '=' in i:
                    (att, value, *_) = i.split('=')
                ann_dict[att] = value

            if feature_type in feature_types:
                table_line = ''
                for attribute in attributes:
                    if attribute in ann_dict:
                        table_line = table_line + ann_dict[attribute] + '\t'
                    else:
                        table_line = table_line + '-' + '\t'
                table_line = table_line + feature_type + '\t' + str(locus_start) + '\t' + str(locus_end) + '\t' + str(strand)
                # Creating snpEff locus id: e.g. CDS_EMRSA15_938198_938683
                # Note: feature.location.start extracted is shifted 1bp leftwards, so +1 need to be added
                snpEff_locus_id = feature_type + '_' + chr_id + '_' + str(locus_start) + '_' + str(locus_end)
                table_line = table_line + '\t' + snpEff_locus_id + '\n'
                print(table_line)
                gff_tab_file.write(table_line)
    in_handle.close()
    gff_tab_file.close()


if __name__ == "__main__":
    _main()
