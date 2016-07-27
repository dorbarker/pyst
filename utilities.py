from Bio import SeqIO
import functools
import os
import tempfile

def basename(path):

    return os.path.splitext(os.path.basename(path))[0]

def format_blast_results(aln):
    
    hsp = aln.hsps[0]
    out = {'AlignmentLength': hsp.aln_length,
           'BitScore': hsp.score,
           'Gaps': hsp.gaps,
           'Mismatches': hsp.positives - hsp.identities,
           'PercentIdentity': hsp.identities / aln.length,
           'QueryAln': hsp.query,
           }

    return out 

def get_query(path):

    with open(path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            return str(record.seq)  # quits after first sequence

def tempdir(func):

    @functools.wraps
    def wrapper(*args, **kwargs):

        with tempfile.TemporaryDirectory() as d:
            return f(*args, **kwargs, temp_dir=d)

    return wrapper
