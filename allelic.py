from Bio.Blast import NCBIXML
from io import StringIO
import subprocess
import utilities
import os
def allelic(alleles, genome):

    def blast_allelic(reference):

        result = blast_reference_allele(reference, genome.path)
        return parse_allelic(result)

    ref = str(alleles.sequences[0].seq)
    exact = find_exact_match(alleles.sequences, genome.sequences)

    return exact or blast_allelic(ref)

def find_exact_match(alleles, genome):

    for contig in genome:

        str_contig = str(contig.seq)  # to allow .index()

        for allele in alleles:

            str_allele = str(allele.seq)  # as above

            try:  # if exact match, return immediately

                res = (str_contig.index(str_allele), str_allele,
                       allele.id, contig.id)
                print(allele.id)
                return res

            except ValueError:
                pass
    else:
        return None  # to be explicit

@utilities.tempdir
def blast_reference_allele(ref, genome_path, temp_dir):

    db_path = os.path.join(temp_dir, utilities.basename(genome_path))
    makeblastdb = ('makeblastdb', '-dbtype', 'nucl',
                   '-in', genome_path, '-out', db_path)

    blastn = ('blastn', '-db', db_path, '-outfmt', '5')

    subprocess.call(makeblastdb, stdout=subprocess.DEVNULL)

    return subprocess.check_output(blastn, input=ref, universal_newlines=True)

def parse_allelic(blast_result):

    result = NCBIXML.read(StringIO(blast_result))
    for aln in result.alignments:
        res = utilities.format_blast_results(aln)
