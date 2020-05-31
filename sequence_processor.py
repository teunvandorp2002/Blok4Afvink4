from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import re


def sequence_type(seq):
    """Uses a regex to determine what type of sequence is given

    :param seq: A string with the input sequence
    :return: The type of sequence and the sequence with just de letters
             (without spaces and enters)
    """
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
    """Transcribes and translates a DNA sequence

    :param seq: A string with the DNA sequence
    :return: A string with the RNA and protein sequences
    """
    dna = str()
    sequence = Seq(seq, IUPAC.unambiguous_dna)
    dna += 'DNA sequence: ' + sequence + '\n'
    dna += 'RNA sequence: ' + str(sequence.transcribe()) + '\n'
    dna += 'Protein sequence: ' + str(sequence.translate(table=1)) + '\n'
    return dna


def blast(seq):
    """Blasts a protein sequence with blastp and writes the results to
    an XML file

    :param seq: A string with the protein sequence
    :return: An XML file with the blast results
    """
    results = NCBIWWW.qblast(program='blastp', database='nr', sequence=seq)
    with open('BLAST_results.xml', "w") as out:
        out.write(results.read())
    results.close()


def parse():
    """Reads the blast XML file and returns the results

    :return: A srting with the results of the blast search
    """
    return_string = str()
    with open("BLAST_results.xml", "r") as blast_file:
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
