from Bio import SeqIO
import functools
import os
import tempfile

def basename(path):

    return os.path.splitext(os.path.basename(path))[0]

def format_blast_results(aln):

    hsp = aln.hsps[0]
    out = {'AlignmentLength': hsp.align_length,
           'BitScore': hsp.score,
           'Gaps': hsp.gaps,
           'Mismatches': hsp.positives - hsp.identities,
           'PercentIdentity': hsp.identities / aln.length,
           'QueryAln': hsp.query,
           }

    return out

def parse_markers(path):

    def convert_typing_test(value):
        lookup = ('pcr', 'allelic', 'repeat',
                  'oligoprobe', 'snp', 'ampliconprobe')

        try:
            return lookup[int(value)]
        except ValueError:
            return value

    d = {}

    with open(path) as f:
        # old style format for now
        first_line = True
        for line in f:

            # skip header
            if first_line:
                first_line = False
                continue

            name, *l = line.strip().split('\t') + ['']  # diag

            d[name] = {'testname': l[0],
                       'testtype': convert_typing_test(l[1]),
                       'forward_primer': l[2],
                       'reverse_primer': l[3],
                       'amplicon_size': int(l[4]),
                       'amplicon_range': float(l[5]),
                       'allelic_db_file': l[6],
                       'repeat_size': l[7]}

    return d

def tempdir(func):

    #@functools.wraps
    def wrapper(*args, **kwargs):

        with tempfile.TemporaryDirectory() as d:
            return func(*args, temp_dir=d, **kwargs)

    return wrapper

class MultiFasta(object):

    def __init__(self, fasta_path):

        self.path = fasta_path
        self.sequences = self.get_sequences()

    def get_sequences(self):

        with open(self.path, 'r') as f:
            return list(SeqIO.parse(f, 'fasta'))
