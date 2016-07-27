#!/usr/bin/env python3

import argparse

def arguments():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-a', '--alleles', help='Path to alleles directory')
    
    parser.add_argument('-t', '--test', required=True,
                        help='Typing test markers file')
    


    return parser.parse_args()

def main():

    args = arguments()

