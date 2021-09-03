"""
How to use it.

> pip install --upgrade wheel
> pip install .
> python demo.py
"""
from pprint import pprint as print

import algorithm

from biosequence.config import AlignmentConfig


def alignment(query, subject, mode=1, return_score=False):
    if mode == 1:
        return algorithm.NeedlemanWunsch(query, subject, return_score,
                                         AlignmentConfig.MATCH,
                                         AlignmentConfig.MISMATCH,
                                         AlignmentConfig.GAP_OPEN,
                                         AlignmentConfig.GAP_EXTEND)
    elif mode == 2:
        return algorithm.SmithWaterman(query, subject, return_score,
                                       AlignmentConfig.MATCH,
                                       AlignmentConfig.MISMATCH,
                                       AlignmentConfig.GAP_OPEN,
                                       AlignmentConfig.GAP_EXTEND)
    else:
        print("Please choose alignment mode:")
        print("1-Global alignment by Needleman-Wunsch")
        print("2-Local alignment by Smith-Waterman")


if __name__ == "__main__":

    import time
    seq1 = "aattttgttttagagacagggtctcagcctcccaagtagctgggacgacaggtgcacaccaccatgcctggctaattttaaaattttttgtagagatggggtcctcacaattttgcccaggctggtcttgaactcctgagctcaagggagcctcctgcctcggcctcccaaagtgttgggattacaggcgtgagcactgcacccagccacctggtctgcgtcttaaaagccttcctgactctcaggactgaaagctgccaccagggcgcctttggaaatcgtcgtaattataacccccccggcctgggcgctgagtccttcccaccagccagcaggcactgaccgggttgcagatcgggagacggaggctcggagaggcccaggggctgctctgccatcccccctttccctgcagcctgggggctccctgacgcctggactccccccctgcaggtgctgcccgagctgcggagatcgtgggcgggcacgaggcgcagccacactcccggccctacatggcctccctgcagatgcgggggaacccgggcagccacttctgcggaggcaccttgatccaccccagcttcgtgctgacggccgcgcactgcctgcgggacatgtgagcggccgcctccacacccctgtccgcccgccccgccctcttcctccagccctggcccggccgctgtccctctgcccggggaggacccagctaagccccgtctgcagaccccaggccccgcgcgcgtgggcagttctggggggaggcccggggcagggtcgccgagggaggggtctggggctgcaccgcggcctcgggaagggccggctgtgggcggcggcgagtgtccagggcgccgaggagtgaccaccccacccccgcagaccccagcgcctggtgaacgtggtgctcggagcccacaacgtgcggacgcaggagcccacccagcagcacttctcggtggctcaggtgtttctgaacaactacgacgcggagaacaaactgaacgacgttctcctcatccaggtgggcgggcagggccgcgagggctcggaggggcacggccagagggctccgggacccccattcctgcagccagcattcattgagcaccactgtatcgcaaccggagcacccactgtataccgggccacgaccgaggtcacgccactgcattccagcctgggtggcagagcaagactccatctcaaaagaaagaaagaaggaaagaaaatgaatgaatacaatagtgacaaatgggacaaagggggtcgtggggcccaggcggagggagcggcatccgcggcgttttgaggtggtgggtgtggtgggtgtggtgggagggcggcccgggcggccaccgtgacctggaagcagcgtctcaccgccgcctgccttctgccccagctgagcagcccagccaacctcagtgcctccgtcgccacagtccagctgccacagcaggaccagccagtgccccacggcacccagtgcctggccatgggctggggccgcgtgggtgcccacgaccccccagcccaggtcctgcaggagctcaatgtcaccgtggtcaccttcttctgccggccacataacatttgcactttcgtccctcgccgcaaggccggcatctgcttcgtaagtaaccgtgcccccaccccgggcaccgggctgccatgaggggaggagggcggcggccagggttccacgcccacctcttagctgtgtggcttcatgctgtgcctcagtctccccacctggatggccgtccctgtcctccagggagactcaggtggccccctgatctgtgatggcatcatccaaggaatagactccttcgtgatctggggatgtgccacccgccttttccctgacttcttcacgcgggtagccctctacgtggactggatccgttccacgctgcgccgtgtggaggccaagggccgcccctgaaccgcccctcccacagcgctggccgggaccccgagcctggctccaaaccctcgaggcggatctttggacagaagcagctcttccccgaacactgtggcgtccgggacggccccacccgtccccccacactccctcccacggggctccgggagacaggccggccctgcacctcaccccaccgtgacctcaataaacgttgaaactcccctggctcctgtctgtccttcctatagg"
    seq2 = "aattttgttttagagacagggtcttgctctgtttgtcctcagcctcccaagtagctgggacgacaggtgcacaccaccatgcctggctaattttaaaattttttgtagagatggggtcctcacaattttgcccaggctggtcttgaactcctgagctcaagggagcctcctgcctcggcctcccaaagtgttgggattacaggcgtgagcactgcacccagccacctggtctgcgtcttaaaagccttcctgactctcaggactgaaagctgccaccagggcgcctttggaaatcgtcgtaattataacccccccggcctgggcgctgagtccttcccaccagccagcaggcactgaccgggttgcagatcgggagacggaggctcggagaggcccaggggctgctctgccatcccccctttccctgcagcctgggggctccctgacgcctggactccccccctgcaggtgctgcccgagctgcggagatcgtgggcgggcacgaggcgcagccacactcccggccctacatggcctccctgcagatgcgggggaacccgggcagccacttctgcggaggcaccttgatccaccccagcttcgtgctgacggccgcgcactgcctgcgggacatgtgagcggccgcctccacacccctgtccgcccgccccgccctcttcctccagccctggcccggccgctgtccctctgcccggggaggacccagctaagccccgtctgcagaccccaggccccgcgcgcgtgggcagttctggggggaggcccggggcagggtcgccgagggaggggtctggggctgcaccgcggcctcgggaagggccggctgtgggcggcggcgagtgtccagggcgccgaggagtgaccaccccacccccgcagaccccagcgcctggtgaacgtggtgctcggagcccacaacgtgcggacgcaggagcccacccagcagcacttctcggtggctcaggtgtttctgaacaactacgacgcggagaacaaactgaacgacgttctcctcatccaggtgggcgggcagggccgcgagggctcggaggggcacggccagagggctccgggacccccattcctgcagccagcattcattgagcaccactgtatcgcaaccggagcacccactgtataccgggccacgaccgaggtcacgccactgcattccagcctgggtggcagagcaagactccatctcaaaagaaagaaagaaggaaagaaaatgaatgaatacaatagtgacaaatgggacaaagggggtcgtggggcccaggcggagggagcggcatccgcggcgttttgaggtggtgggtgtggtgggtgtggtgggagggcggcccgggcggccaccgtgacctggaagcagcgtctcaccgccgcctgccttctgccccagctgagcagcccagccaacctcagtgcctccgtcgccacagtccagctgccacagcaggaccagccagtgccccacggcacccagtgcctggccatgggctggggccgcgtgggtgcccacgaccccccagcccaggtcctgcaggagctcaatgtcaccgtggtcaccttcttctgccggccacataacatttgcactttcgtccctcgccgcaaggccggcatctgcttcgtaagtaaccgtgcccccaccccgggcaccgggctgccatgaggggaggagggcggcggccagggttccacgcccacctcttagctgtgtggcttcatgctgtgcctcagtctccccacctggatggccgtccctgtcctccagggagactcaggtggccccctgatctgtgatggcatcatccaaggaatagactccttcgtgatctggggatgtgccacccgccttttccctgacttcttcacgcgggtagccctctacgtggactggatccgttccacgctgcgccgtgtggaggccaagggccgcccctgaaccgcccctcccacagcgctggccgggaccccgagcctggctccaaaccctcgaggcggatctttggacagaagcagctcttccccgaacactgtggcgtccgggacggccccacccgtccccccacactccctcccacggggctccgggagacaggccggccctgcacctcaccccaccgtgacctcaataaacgttgaaactcccctggctcctgtctgtccttcctatagg"

    print(f"序列长度{len(seq1)} x {len(seq2)}")
    start_time = time.time()
    result = alignment(seq1, seq2)
    print(result)
    end_time = time.time()
    print("C: 总共用时%s" % (end_time - start_time))
    start_time = time.time()
    result = alignment(seq1, seq2, mode=2)
    print(result)
    end_time = time.time()
    print("Python: 总共用时%s" % (end_time - start_time))
