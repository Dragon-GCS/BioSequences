test_dna = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'


################################## Coden Table ##################################
TABLE = {
        "AAA":"K",	"AAC":"N",	"AAG":"K",	"AAU":"N",	"ACA":"T",	"ACC":"T",	"ACG":"T",	"ACU":"T",	
        "AGA":"R",	"AGC":"S",	"AGG":"R",	"AGU":"S",	"AUA":"I",	"AUC":"I",	"AUG":"M",	"AUU":"I",	
        "CAA":"Q",	"CAC":"H",	"CAG":"Q",	"CAU":"H",	"CCA":"P",	"CCC":"P",	"CCG":"P",	"CCU":"P",	
        "CGA":"R",	"CGC":"R",	"CGG":"R",	"CGU":"R",	"CUA":"L",	"CUC":"L",	"CUG":"L",	"CUU":"L",	
        "GAA":"E",	"GAC":"D",	"GAG":"E",	"GAU":"D",	"GCA":"A",	"GCC":"A",	"GCG":"A",	"GCU":"A",	
        "GGA":"G",	"GGC":"G",	"GGG":"G",	"GGU":"G",	"GUA":"V",	"GUC":"V",	"GUG":"V",	"GUU":"V",	
        "UAA":"Stop",	"UAC":"Y",	"UAG":"Stop",	"UAU":"Y",	"UCA":"S",	"UCC":"S",	"UCG":"S",	"UCU":"S",	
        "UGA":"Stop",	"UGC":"C",	"UGG":"W",	"UGU":"C",	"UUA":"L",	"UUC":"F",	"UUG":"L",	"UUU":"F",
        }

START_CODON = "AUG"

def setTable(json_file):
    import json
    from pathlib import Path

    global TABLE

    filename = Path(__file__).resolve().parent.joinpath("coden_table", json_file)
    with open(filename) as f:
        TABLE = json.load(f)

################################ Align Paremeter ################################
class AlignmentConfig:
    MATCH = 2
    MISMATCH = -3
    GAP_OPEN = -3
    GAP_EXTEND = -3

def setAlignPara(match=2, mismatch=-3, gap_open=-3, gap_extend=-3):
    AlignmentConfig.MATCH = match
    AlignmentConfig.MISMATCH = mismatch
    AlignmentConfig.GAP_OPEN = gap_open
    AlignmentConfig.GAP_EXTEND = gap_extend

############################### Molecular Weight ################################
# All molecular weight were subtracted a H2O（18.0）
MW = {
    "Peptide_MW":{
            "A":71.1,   "C":103.2,  "D":115.1,   "E":129.1,  "F":147.2, 	
            "G":57.1,   "H":137.2,  "I":113.2,   "K":128.2,  "L":113.2, 	
            "M":131.2,  "N":114.1,  "P":97.1,    "Q":128.15, "R":156.20, 	
            "S":87.09,  "T":101.16, "V":99.15,   "W":186.22, "Y":163.19,
            },
    "DNA_MW": {"A":313.2,  "C":289.2,  "G":329.2,  "T":304.2,},
    "RNA_MW": {"A":329.2,  "C":305.2,  "G":345.2,  "U":306.2,},
    
}
# refer to https://www.thermofisher.cn/cn/zh/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html

############################## Nuclear Acid Info ################################
NC_INFO = {
    "DNA_COMPLEMENT": {"A":"T", "C":"G", "G":"C", "T":"A",},
    "RNA_COMPLEMENT": {"A":"T", "C":"G", "G":"C", "U":"A",}
    }


############################## Peptide Info ################################
# Author(s): Kyte J., Doolittle R.F.
# Reference: J. Mol. Biol. 157:105-132(1982).
HYDROPATHY = {
                "A": 1.8,   "C": 2.5,   "D": -3.5,  "E": -3.5,  "F": 2.8, 
                "G": -0.4,  "H": -3.2,  "I": 4.5,   "K": -3.9,  "L": 3.8, 
                "M": 1.9,   "N": -3.5,  "P": -1.6,  "Q": -3.5,  "R": -4.5,
                "S": -0.8,  "T": -0.7,  "V": 4.2,   "W": -0.9,  "Y": -1.3, 
                }
# value from "EMBOSS"
pK = {
        "K":10.8,   "R":12.5,  "H":6.5,  "D":3.9,  "E":4.1,
        "C":8.5,    "Y":10.1,   "N_term":8.6,   "C_term": 3.6,
    }