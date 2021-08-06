# coding = utf-8
# Dragon's Python3.8 code
# Created at 2021/08/06 16:19:00
# Edit with VS Code

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


def get_next(string):
    """
    Find the next-array of pattern string.
    Args:
        string: the pattern string
    Returns:
        next: the next-array of pattern string
    """
    length = len(string)
    next = [0]*length
    next[0] = -1
    i = 0
    index = -1

    while i < length-1:
        if index==-1 or string[index] == string[i]:
            i += 1
            index += 1
            if string[index] == string[i]:
                next[i] = next[index]
            else:
                next[i] = index
        else:
            index = next[index]

    return next


def kmp_search(s1, s2):
    """
    Search s2 in s1 and ruturn index if s2 in s1.
    Args:
        s1: main string be searched
        s2: pattern string
    Returns:
        the index of s2 pair to s1
    """
    next = get_next(s2)
    i = 0
    j = 0

    while i < len(s1) and j < len(s2):
        if j==-1 or s1[i] == s2[j]:
            i += 1
            j += 1
        else:
            j = next[j]

    if j == len(s2):
        return i - j

    return None


if __name__ == '__main__':
    def read_fasta_test():
        seqs, ids = read_fasta("test.txt")
        assert seqs != []
        assert seqs[0] != ""
        assert len(seqs) == len(ids)

    def kmp_test():
        assert kmp_search("ababad", 'abad') == 2
    
    read_fasta_test()
    kmp_test()