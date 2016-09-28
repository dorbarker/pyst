#!/usr/bin/env python3

import argparse
import utilities
import allelic
import os
import time

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alleles', help='Path to alleles directory')

    parser.add_argument('-t', '--test', required=True,
                        help='Typing test markers file')

    parser.add_argument('genome', help='FASTA formatted genome')

    return parser.parse_args()

@utilities.tempdir
def run_test(markers, allele_dir, genome_path, temp_dir):

    genome = utilities.MultiFasta(genome_path)

    for m in markers:
        if markers[m]['testtype'] == 'allelic':

            allele_path = os.path.join(allele_dir, markers[m]['allelic_db_file'])
            alleles = utilities.MultiFasta(allele_path)

            res = allelic.allelic(alleles, genome, temp_dir)
            print(allelic.match_allele(res['SubjAln'], alleles.str_seqs))
def main():

    args = arguments()
    t1 = time.time()
    markers = utilities.parse_markers(args.test)
    run_test(markers, args.alleles, args.genome)

    print(time.time() - t1)

if __name__ == '__main__':
    main()
