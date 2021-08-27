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
