TABLE = {
        "AAA":"K",	"AAC":"N",	"AAG":"K",	"AAT":"N",	"ACA":"T",	"ACC":"T",	"ACG":"T",	"ACT":"T",	
        "AGA":"R",	"AGC":"S",	"AGG":"R",	"AGT":"S",	"ATA":"I",	"ATC":"I",	"ATG":"M",	"ATT":"I",	
        "CAA":"Q",	"CAC":"H",	"CAG":"Q",	"CAT":"H",	"CCA":"P",	"CCC":"P",	"CCG":"P",	"CCT":"P",	
        "CGA":"R",	"CGC":"R",	"CGG":"R",	"CGT":"R",	"CTA":"L",	"CTC":"L",	"CTG":"L",	"CTT":"L",	
        "GAA":"E",	"GAC":"D",	"GAG":"E",	"GAT":"D",	"GCA":"A",	"GCC":"A",	"GCG":"A",	"GCT":"A",	
        "GGA":"G",	"GGC":"G",	"GGG":"G",	"GGT":"G",	"GTA":"V",	"GTC":"V",	"GTG":"V",	"GTT":"V",	
        "TAA":"Stop",	"TAC":"Y",	"TAG":"Stop",	"TAT":"Y",	"TCA":"S",	"TCC":"S",	"TCG":"S",	"TCT":"S",	
        "TGA":"Stop",	"TGC":"C",	"TGG":"W",	"TGT":"C",	"TTA":"L",	"TTC":"F",	"TTG":"L",	"TTT":"F",
        }
def setTable(json_file):
    import json
    from pathlib import Path

    global TABLE

    filename = Path(__file__).resolve().parent.joinpath("coden_table", json_file)
    with open(filename) as f:
        TABLE = json.load(f)


START_CODON = "ATG"
 
test_dna = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'

class AlignmentConfig:
    MATCH = 2
    MISMATCH = -3
    GAP_OPEN = -3
    GAP_EXTEND = -3

def setAlignPara(match, mismatch, gap_open, gap_extend):
    global AlignmentConfig
    AlignmentConfig.MATCH = match
    AlignmentConfig.MISMATCH = mismatch
    AlignmentConfig.GAP_OPEN = gap_open
    AlignmentConfig.GAP_EXTEND = gap_extend