from typing import Any, Dict, List, Optional, Union

############################### About utils ##################################
#: :meta hide-value: SYMBOL["printAlign"] - Use to print alignment sequences
SYMBOL: Dict[str, List[str]] = {
    # replace "| · -" or other fixed width character
    "printAlign": ["┃", "•", "━"]
}

################################## Codon Table ##################################
#: :meta hide-value: mRNA codon table, used to ``transcript()``
CODON_TABLE: Dict[str, str] = {
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
#: :meta hide-value: | Start code list, default is ["AUG"], used to ``get_orf()``.
#: | All codon in this list will be set to start of ORF.
#: | When ``from bioseq.utils import START_CODON``, to change a start codon, don't assign to ``START_CODON``, use ``START_CODON[0] = "ATT"`` to instead.
START_CODON: List[str] = ["AUG"]

################################ Align Parameter ################################


class AlignmentConfig:
    """Align parameters

    Attributes:
        MATCH (float): score when meet a match pair
        MISMATCH (float): score when meet a mismatch pair
        GAP_OPEN (float): score when open a gap
        GAP_EXTEND (float): score when extend a gap
    """
    MATCH: float = 2
    MISMATCH: float = -3
    GAP_OPEN: float = -3
    GAP_EXTEND: float = -3

############################### Molecular Weight ################################


#: :meta hide-value: | Molecular weight table, include "Peptide", "DNA", "RNA".
#: | reference https://www.thermofisher.cn/cn/zh/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
MW: Dict = {
    "Peptide": {
        "A": 89.1,  "C": 121.2, "D": 133.1, "E": 147.1, "F": 165.2,
        "G": 75.1,  "H": 155.2, "I": 131.2, "K": 146.2, "L": 131.2,
        "M": 149.2, "N": 132.1, "P": 115.1, "Q": 146.2, "R": 174.2,
        "S": 105.1, "T": 119.1, "V": 117.1, "W": 204.2, "Y": 181.2,
    },
    "DNA": {"A": 313.2, "C": 289.2, "G": 329.2, "T": 304.2, },
    "RNA": {"A": 329.2, "C": 305.2, "G": 345.2, "U": 306.2, },
}

############################## Nuclear Acid Info ################################
#: :meta hide-value: Complementary table of nucleic acid, include "DNA_COMPLEMENT", "RNA_COMPLEMENT".
NC_INFO: Dict = {
    "DNA_COMPLEMENT": {"A": "T", "C": "G", "G": "C", "T": "A", },
    "RNA_COMPLEMENT": {"A": "U", "C": "G", "G": "C", "U": "A", }
}

################################# Peptide Info ##################################
#: :meta hide-value: | Hydrophobicity of 20 amino acids, the higher, the more hydrophobic
#: | Author(s): Kyte J., Doolittle R.F.
#: | Reference: J. Mol. Biol. 157:105-132(1982).
HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
#: :meta hide-value: | pK value of charged amino acids.
#: | reference ``biopython.SeqUtils.IsoelectricPoint``, Value from EMBOSS Database
PK = {
    "Nterm": 8.6, "Cterm": 3.6,
    "pos_pK": {"K": 10.8, "R": 12.5, "H": 6.5},
    "neg_pK": {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
}
