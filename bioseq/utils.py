from typing import Iterable, List, Union
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.error import HTTPError

from bioseq.config import SYMBOL
from bioseq import DNA, RNA, Peptide, Sequence

EUTILS_POST = {
    "db": "",           # 数据库
    "rettype": "fasta", # 数据类型
    "retmode": "text",  # 返回类型
    "id": "",           # uid
}
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"


def fetchNCBI(uid: str) -> Union[DNA, RNA, Peptide]:
    """
    Fetch sequence corresponding to UID from NCBI E-utilities. Only support RNA, mRNA(DNA), Protein.

    NM_(mRNA):      Protein-coding transcripts (usually curated)
    NR_(RNA ):      Non-protein-coding transcripts
    XM_(mRNA):      Predicted model protein-coding transcript
    XR_(RNA ):      Predicted model non-protein-coding transcript
    AP_(Protein):   Annotated on AC_ alternate assembly
    NP_(Protein):   Associated with an NM_ or NC_ accession
    YP_(Protein):   Annotated on genomic molecules without an instantiated transcript record
    XP_(Protein):   Predicted model, associated with an XM_ accession
    WP_(Protein):   Non-redundant across multiple strains and species

    Args:
        uid: NCBI's unique id
    Returns:
        A sequence object corresponding to UID
    """
    # check uid
    if uid[:3] in ["AP_", "NP_", "YP_", "XP_", "WP_"]:
        sequence = Peptide()
        EUTILS_POST["db"] = "protein"
    elif uid[:3] in ["NM_", "XM_"]:
        EUTILS_POST["db"] = "nuccore"
        sequence = DNA()
    elif uid[:3] in ["NR_", "XR_"]:
        EUTILS_POST["db"] = "nuccore"
        sequence = RNA()
    else:
        raise ValueError(f"{uid} is not a support uid")

    EUTILS_POST["id"] = uid
    try:
        raw_info: List[str] = urlopen(
            EUTILS_URL + urlencode(EUTILS_POST)).read().decode().split("\n")
        sequence._seq, sequence.info = "".join(
            raw_info[1:]), raw_info[0].lstrip(">")
    except HTTPError as e:
        print(e)
    finally:
        return sequence


def loadFasta(filename: str) -> Iterable[Sequence]:
    """
    Read fasta file

    Args:
        filename: the fasta file's name
    Returns:
        seq_list: the list of sequences
        seq_ids: the list of sequence's ids
    """
    seq = ""
    info = ""

    with open(filename) as f:
        while line := f.readline():
            # process the last line
            if not line.strip():
                yield Sequence(seq, info)

            if line.startswith(">"):
                # yield previous sequence
                if seq:
                    yield Sequence(seq, info)
                    seq = ""    # reset seq
                info = line[1:].strip()
            else:
                seq += line.strip()

        yield Sequence(seq, info)


def printAlign(
        sequence1: str,
        sequence2: str,
        spacing: int = 10,
        line_width: int = 30,
        show_seq: bool = True) -> None:
    """
    Print two sequence by a pretty format

    Args:
        sequence1:  Sequence1
        sequence2:  Sequence2
        spacing:    A space each $spacing char
        line_width: the width of each line
        show_seq:   if False, only print the alignment result
    """
    symbol_line = ""
    format_seq1 = ""
    format_seq2 = ""
    count = 0
    length = len(sequence1) if len(sequence1) < len(
        sequence2) else len(sequence2)
    match_symbol, mismatch_symbol, gap_symbol = SYMBOL["printAlign"]

    for i in range(length):
        base1 = sequence1[i]
        base2 = sequence2[i]
        if base1 == base2:
            symbol_line += match_symbol
        elif base1 == "-" or base2 == "-":
            symbol_line += gap_symbol
        else:
            symbol_line += mismatch_symbol

        format_seq1 += base1
        format_seq2 += base2
        count += 1

        if count % spacing == 0:
            format_seq1 += " "
            format_seq2 += " "
            symbol_line += " "
        if count % line_width == 0:
            count = 0

    space_num = 0
    space_each_line = line_width // spacing

    for i in range(0, length, line_width):
        start = i + space_num
        end = start + line_width + space_each_line
        space_num += space_each_line
        if show_seq:
            print(f"{i + 1:>5}", end=" ")
            print(format_seq1[start: end])

            print(f" " * 5, end=" ")
            print(symbol_line[start: end])

            print(f"{i + 1:>5}", end=" ")
            print(format_seq2[start: end])
        else:
            print(f"{i + 1:>5}", end=" ")
            print(symbol_line[start: end])

        print()
