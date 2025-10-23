#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio import SeqFeature as sf
# from BCBio import GFF

# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to define loci in a reference genome that can be used as units of association in a GWAS.
# CDS are extracted from an input GFF and added a chosen downstream region (putative regulatory region).
# NOTE: define_loci_in_reference_genome.manual.py script has the same funcionality as define_loci_in_reference_genome.py
# but does not make use of BCBio to parse the input GFF. Instead, define_loci_in_reference_genome.manual.py had to be
# written to deal with the format of Mycobacterium_tuberculosis_H37Rv_gff_v4.gff

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to define loci in a reference genome"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-i", "--input_fasta", action="store", dest="input_fasta",
        help="Reference genome FASTA file",
        required=True, metavar="INPUT_FASTA"
    )
    group.add_argument(
        "-g", "--input_gff", action="store", dest="input_gff",
        help="Reference genome GFF file",
        required=True, metavar="INPUT_GFF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table with defined loci in the reference genome",
        required=True, metavar="OUTPUT_TABLE"
    )
    group = parser.add_argument_group('loci defining arguments')
    group.add_argument(
        "-d", "--downstream_bp", action="store", dest="downstream_bp",
        help="Downstream base pairs from CDS to be considered as regulatory region (e.g. 200)",
        required=True, metavar="DNSTRM_BP"
    )
    group.add_argument(
        "-l", "--gff_locus_id_field", action="store", dest="gff_locus_id_field",
        help="Field in the GFF file that contain locus Ids (e.g. locus_tag or gene)",
        required=True, metavar="LOCUS_ID"
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
    input_files = [args.input_fasta, args.input_gff]
    for input_file in input_files:
        if not os.path.isfile(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Reading input FASTA file
    logging.info(f"Reading sequence's Ids and lengths from {args.input_fasta}")
    fasta_records_length = dict()
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        fasta_records_length[record.id] = len(record)
        print(record.id + ' ' + str(len(record)))

    # Defining feature types and attributes to extract
    feature_types = ['CDS']
    # feature_types = ['gene']
    gff_locus_id_field = args.gff_locus_id_field

    # Dictionary to save locus annotation per chromosome position
    positions_annotated = dict()

    # Reading locus annotation from GFF
    logging.info(f"Opening input annotated genome {args.input_gff}")
    in_handle = open(args.input_gff)
    for line in in_handle:
        # print(line)
        items = line.strip().split('\t')
        feature_type = items[2]
        if feature_type in feature_types:
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
            # print('ann_dict[Locus] ' + ann_dict['Locus'])
            if gff_locus_id_field not in ann_dict:
                logging.error(f"{gff_locus_id_field} not found in GFF line: {line}")
                exit(0)
            locus_id = ann_dict[gff_locus_id_field]
            for pos in range(int(locus_start), int(locus_end)):
                if positions_annotated.get(pos) is None:
                    positions_annotated[pos] = {}
                if positions_annotated[pos].get("locus_start") is None:
                    positions_annotated[pos]["locus_start"] = str(locus_start)
                else:
                    positions_annotated[pos]["locus_start"] = positions_annotated[pos]["locus_start"] + ';' + str(locus_start)
                if positions_annotated[pos].get("locus_end") is None:
                    positions_annotated[pos]["locus_end"] = str(locus_end)
                else:
                    positions_annotated[pos]["locus_end"] = positions_annotated[pos]["locus_end"] + ';' + str(locus_end)
                if positions_annotated[pos].get("locus_strand") is None:
                    positions_annotated[pos]["locus_strand"] = str(strand)
                else:
                    positions_annotated[pos]["locus_strand"] = positions_annotated[pos]["locus_strand"] + ';' + str(strand)
                if positions_annotated[pos].get("locus_id") is None:
                    positions_annotated[pos]["locus_id"] = str(locus_id) + ';'
                else:
                    positions_annotated[pos]["locus_id"] = positions_annotated[pos]["locus_id"] + str(locus_id) + ';'
    in_handle.close()

    # Creating table for all chromosome positions
    header = 'chromosome\tposition\tlocus_id\tlocus_start\tlocus_end\tlocus_strand\n'
    output_table_file = open(args.output_table, 'w')
    output_table_file.write(header)

    for chromosome in fasta_records_length:
        for position in range(1, fasta_records_length[chromosome]):
            output_line = chromosome + '\t' + str(position)
            if position in positions_annotated:
                output_line = output_line + '\t' + positions_annotated[position]['locus_id'] + '\t'\
                              + positions_annotated[position]['locus_start'] + '\t' + positions_annotated[position]["locus_end"] + '\t'\
                              + positions_annotated[position]["locus_strand"]
            else:
                output_line = output_line + '\t' + '\t'.join("-" * 4)
            print(output_line)
            output_table_file.write(output_line + '\n')
    output_table_file.close()


if __name__ == "__main__":
    _main()


