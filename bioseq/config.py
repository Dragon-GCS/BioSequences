class TestSeq:
    name = "Homo sapiens actin beta (ACTB), mRNA"
    ID = "NM_001101.5"
    DNA = "atggatgatgatatcgccgcgctcgtcgtcgacaacggctccggcatgtgcaaggccggcttcgcgggcg" \
          "acgatgccccccgggccgtcttcccctccatcgtggggcgccccaggcaccagggcgtgatggtgggcat" \
          "gggtcagaaggattcctatgtgggcgacgaggcccagagcaagagaggcatcctcaccctgaagtacccc" \
          "atcgagcacggcatcgtcaccaactgggacgacatggagaaaatctggcaccacaccttctacaatgagc" \
          "tgcgtgtggctcccgaggagcaccccgtgctgctgaccgaggcccccctgaaccccaaggccaaccgcga" \
          "gaagatgacccagatcatgtttgagaccttcaacaccccagccatgtacgttgctatccaggctgtgcta" \
          "tccctgtacgcctctggccgtaccactggcatcgtgatggactccggtgacggggtcacccacactgtgc" \
          "ccatctacgaggggtatgccctcccccatgccatcctgcgtctggacctggctggccgggacctgactga" \
          "ctacctcatgaagatcctcaccgagcgcggctacagcttcaccaccacggccgagcgggaaatcgtgcgt" \
          "gacattaaggagaagctgtgctacgtcgccctggacttcgagcaagagatggccacggctgcttccagct" \
          "cctccctggagaagagctacgagctgcctgacggccaggtcatcaccattggcaatgagcggttccgctg" \
          "ccctgaggcactcttccagccttccttcctgggcatggagtcctgtggcatccacgaaactaccttcaac" \
          "tccatcatgaagtgtgacgtggacatccgcaaagacctgtacgccaacacagtgctgtctggcggcacca" \
          "ccatgtaccctggcattgccgacaggatgcagaaggagatcactgccctggcacccagcacaatgaagat" \
          "caagatcattgctcctcctgagcgcaagtactccgtgtggatcggcggctccatcctggcctcgctgtcc" \
          "accttccagcagatgtggatcagcaagcaggagtatgacgagtccggcccctccatcgtccaccgcaaat" \
          "gcttctag"
    Peptide = "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTL" \
              "KYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQ" \
              "AVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAERE" \
              "IVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHET" \
              "TFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILA" \
              "SLSTFQQMWISKQEYDESGPSIVHRKCF"


############################### About utils ##################################

SYMBOL = {
        "printAlign": ("┃", "•", "━")  # replace "| · -" or other fixed width character
}

################################## Coden Table ##################################
TABLE = {
        "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
        "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
        "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
        "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
        "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
        "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
        "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
        "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y",
        "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
        "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C",
        "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F"
}

START_CODON = ["AUG"]


def setStartCoden(coden):
    global START_CODON
    if isinstance(coden, str) and len(coden) == 3:
        START_CODON = [coden]
    elif isinstance(coden, list):
        START_CODON = coden
    else:
        raise ValueError("Coden should be a str or list of str")


################################ Align Parameter ################################
class AlignmentConfig:
    MATCH:float = 2
    MISMATCH:float = -3
    GAP_OPEN:float = -3
    GAP_EXTEND:float = -3


def setAlignPara(match:float = 2, mismatch:float = -3, gap_open:float = -3, gap_extend:float = -3):
    AlignmentConfig.MATCH = match
    AlignmentConfig.MISMATCH = mismatch
    AlignmentConfig.GAP_OPEN = gap_open
    AlignmentConfig.GAP_EXTEND = gap_extend


############################### Molecular Weight ################################
# reference https://www.thermofisher.cn/cn/zh/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and
# -rna-molecular-weights-and-conversions.html
# All molecular weight were subtracted a H2O（18.0）
MW = {
        "Peptide_MW": {
                "A": 71.1, "C": 103.2, "D": 115.1, "E": 129.1, "F": 147.2,
                "G": 57.1, "H": 137.2, "I": 113.2, "K": 128.2, "L": 113.2,
                "M": 131.2, "N": 114.1, "P": 97.1, "Q": 128.15, "R": 156.20,
                "S": 87.09, "T": 101.16, "V": 99.15, "W": 186.22, "Y": 163.19,
        },
        "DNA_MW"    : {"A": 313.2, "C": 289.2, "G": 329.2, "T": 304.2, },
        "RNA_MW"    : {"A": 329.2, "C": 305.2, "G": 345.2, "U": 306.2, },
}

############################## Nuclear Acid Info ################################
NC_INFO = {
        "DNA_COMPLEMENT": {"A": "T", "C": "G", "G": "C", "T": "A", },
        "RNA_COMPLEMENT": {"A": "U", "C": "G", "G": "C", "U": "A", }
}

################################# Peptide Info ##################################
# Author(s): Kyte J., Doolittle R.F.
# Reference: J. Mol. Biol. 157:105-132(1982).
HYDROPATHY = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
# reference biopython.SeqUtils.IsoelectricPoint
# data from EMBOSS
PK = {
        "Nterm" : 8.6, "Cterm": 3.6,
        "pos_pK": {"K": 10.8, "R": 12.5, "H": 6.5},
        "neg_pK": {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
}
