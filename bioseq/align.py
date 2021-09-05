from bioseq import algorithm
from bioseq.config import AlignmentConfig


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
