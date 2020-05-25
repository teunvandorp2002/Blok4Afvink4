from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import re


def sequence_type(seq):
    seq = ''.join(filter(str.isalpha, seq)).upper()
    if re.search("^[ACTG]*$", seq):
        return "DNA", seq
    elif re.search("^[AUCG]*$", seq):
        return "RNA", seq
    elif re.search("^[^BJXZ]*$", seq):
        return "Eiwit", seq
    else:
        return "Geen DNA, RNA of eiwit", seq


def dna_processor(seq):
    dna = str()
    sequence = Seq(seq, IUPAC.unambiguous_dna)
    dna += 'DNA sequence: ' + sequence + '\n'
    dna += 'RNA sequence: ' + str(sequence.transcribe()) + '\n'
    dna += 'Protein sequence: ' + str(sequence.translate(table=1)) + '\n'
    return dna


def blast(seq):
    results = NCBIWWW.qblast(program='blastp', database='nr', sequence=seq)
    with open('BLAST_results.bxml', "w") as out:
        out.write(results.read())
    results.close()


def parse():
    return_string = str()
    with open("BLAST_results.bxml", "r") as blast_file:
        result = NCBIXML.read(blast_file)
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                return_string += '-' * 100
                return_string += '\nHit: ' + str(alignment.title)
                return_string += '\nLength: ' + str(alignment.length)
                return_string += '\ne value: ' + str(hsp.expect) + '\n'
    if return_string == "":
        return "Geen matches gevonden"
    return return_string
