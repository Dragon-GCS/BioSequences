from bioseq.config import SYMBOL


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


def printAlign(sequence1, sequence2, spacing = 10, line_width = 30, show_seq = True):
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
    length = len(sequence1) if len(sequence1) < len(sequence2) else len(sequence2)
    match_symbol, mismatch_symbol, gap_symbol = SYMBOL["printAlign"]

    for i in range(length):
        base1 = sequence1[i]
        base2 = sequence2[i]
        if base1 == base2:
            symbol_line += match_symbol
        elif base1 == "-" or  base2 == "-":
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
            print(format_seq1[start : end])
        
            print(f" " * 5, end=" ")
            print(symbol_line[start : end])

            print(f"{i + 1:>5}", end=" ")
            print(format_seq2[start : end])
        else:
            print(f"{i + 1:>5}", end=" ")
            print(symbol_line[start : end])
        
        print()