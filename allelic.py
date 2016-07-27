from Bio.Blast import NCBIXML
from io import StringIO
import subprocess
import utilities

def find_exact_match(alleles, genome):

    for contig in genome:

        str_contig = str(contig.seq)  # to allow .index()

        for allele in alleles:

            str_allele = str(contig.seq)  # as above

            try:  # if exact match, return immediately 
            
                res = (str_contig.index(str_allele), str_allele,
                       allele.id, contig.id)
                
                return res 
            
            except IndexError:
                pass
    else:
        return None  # to be explicit

@utilities.tempdir
def blast_reference_allele(alleles_path, genome_path, temp_dir):

    db_path = os.path.join(temp_dir, utilities.basename(genome_path))

    makeblastdb = ('makeblastdb', '-dbtype', 'nucl',
                   '-in', genome_path, '-out', db_path)

    blastn = ('blastn', '-db', db_path, '-outfmt', '5')

    query = utilities.get_query(alleles_path)

    subprocess.call(makeblastdb)

    return subprocess.check_output(blastn, input=query, universal_newlines=True)

def parse_allelic(blast_result):
    
    for aln in NCBIXML(StringIO(blast_result)).alignments: 
        res = utilities.format_blast_results(aln)  
