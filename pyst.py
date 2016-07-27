#!/usr/bin/env python3

import argparse
import utilities
import allelic
import os

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alleles', help='Path to alleles directory')

    parser.add_argument('-t', '--test', required=True,
                        help='Typing test markers file')

    parser.add_argument('genome', help='FASTA formatted genome')

    return parser.parse_args()

def run_test(markers, allele_dir, genome_path):

    genome = utilities.MultiFasta(genome_path)

    for m in markers:
        if markers[m]['testtype'] == 'allelic':

            allele_path = os.path.join(allele_dir, markers[m]['allelic_db_file'])
            alleles = utilities.MultiFasta(allele_path)

            res = allelic.allelic(alleles, genome)

def main():

    args = arguments()
    markers = utilities.parse_markers(args.test)
    run_test(markers, args.alleles, args.genome)

if __name__ == '__main__':
    main()
