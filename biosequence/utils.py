from typing import Tuple, Union, List
from biosequence.sequence import Sequence
from biosequence.config import SYMBOL

def read_fasta(filename):
    """
    Read content of fasta file.
    Args:
        filename: the fasta file's name
    Returns:
        seq_list: the list of sequences
        seq_ids: the list of sequence's ids
    """
    seq = ""
    seq_list =[]
    seq_ids = []

    with open(filename) as f:
        for line in f.readlines():
            
            if line.startswith(">"):
                seq_ids.append(line[1:].strip())

                if seq:
                    seq_list.append(seq)
                    seq = ""
                continue

            seq += line.strip()

        seq_list.append(seq)

    return seq_list, seq_ids 

def printAlign(sequence1: Union[str, Sequence], sequence2:Union[str, Sequence],spacing:int = 10, line_width:int = 30, show_sequence:bool = True) -> None:
    """
    Print two sequence by a pretty format
    Args:
        sequence1
        sequence2
        line_width: the num of each line
    """
    symbol_line = ""
    format_seq1 = ""
    format_seq2 = ""
    length = len(sequence1) if len(sequence1) < len(sequence2) else len(sequence2)
    match_symbol, mismathc_symbol, gap_symbol = SYMBOL["printAlign"].values()

    for i in range(0, length):
        if (base1 := sequence1[i]) == (base2 := sequence2[i]):
            symbol_line += match_symbol
        elif base1 == "-" or  base2 == "-":
            symbol_line += gap_symbol
        else:
            symbol_line += mismathc_symbol


        format_seq1 += base1
        format_seq2 += base2

        if (i + 1) % spacing == 0:
            format_seq1 += " "
            format_seq2 += " "
            symbol_line += " "

    space_num = 0
    space_each_line = line_width // spacing

    for i in range(0, length, line_width):
        start = i + space_num
        end = start + line_width + space_each_line
        space_num += space_each_line
        if show_sequence:
            print(f"{i + 1:>5}", end=" ")
            print(format_seq1[start : end])
        
            print(f" " * 5, end=" ")
            print(symbol_line[start : end])

            print(f"{i + 1:>5}", end=" ")
            print(format_seq2[start : end])
        else:
            print(f"{i + 1:>5}", end=" ")
            print(symbol_line[start : end])
        
        print()
