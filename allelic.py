from Bio.Blast import NCBIXML
from io import StringIO
import subprocess
import utilities
import os

def allelic(alleles, genome,  temp_dir):

    ref = str(alleles.sequences[0].seq)
    result = blast_reference_allele(ref, genome.path, temp_dir)

    return parse_allelic(result)

def blast_reference_allele(ref, genome_path, temp_dir):

    db_path = os.path.join(temp_dir, utilities.basename(genome_path))

    makeblastdb = ('makeblastdb', '-dbtype', 'nucl',
                   '-in', genome_path, '-out', db_path)

    blastn = ('blastn', '-db', db_path, '-word_size', '7', '-outfmt', '5')

    try:
        return subprocess.check_output(blastn, input=ref, universal_newlines=True)

    except:
        subprocess.call(makeblastdb, stdout=subprocess.DEVNULL)
        return subprocess.check_output(blastn, input=ref, universal_newlines=True)

def parse_allelic(blast_result):

    result = NCBIXML.read(StringIO(blast_result))
    for aln in result.alignments:
        yield utilities.format_blast_results(aln)

def match_allele(match, alleles):

    try:
        return alleles.index(match.replace('-', '')) + 1
    except ValueError:
        return -1
