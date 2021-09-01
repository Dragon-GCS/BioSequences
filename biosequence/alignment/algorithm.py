import platform     # detect system for using .dll or .so to boost alignment
import os
from functools import wraps
from ctypes import *
from biosequence.config import AlignmentConfig

def cAlgorithm(func):
    @wraps(func)
    def comprise(query, subject, return_score=False):
        if not isinstance(query, bytes):
            query = bytes(query, encoding="utf-8")

        if not isinstance(subject, bytes):
            subject = bytes(subject, encoding="utf-8")

        query = create_string_buffer(query)
        subject = create_string_buffer(subject)
        aligned_query = create_string_buffer(b"", len(query) + len(subject))
        aligned_subject = create_string_buffer(b"", len(query) + len(subject))
        score = c_float()
        match = c_float(AlignmentConfig.MATCH)
        mismatch = c_float(AlignmentConfig.MISMATCH)
        gap_open = c_float(AlignmentConfig.GAP_OPEN)
        gap_extend = c_float(AlignmentConfig.GAP_EXTEND)

        if platform.system() == "Windows":
            suffix = ".dll"
        elif platform.system() == "Linux":
            suffix = ".so"
        else:
            suffix = ".dll"
            print("Can't detect system, using .dll to boost Alignment")

        dll_path = os.path.join(os.path.dirname(__file__), "algorithm" + suffix)

        cAlign =func(dll_path)
        cAlign(query, subject, aligned_query, aligned_subject, pointer(score), match, mismatch, gap_open, gap_extend)

        query = str(aligned_query.value, encoding="utf-8")
        subject = str(aligned_subject.value, encoding="utf-8")

        if return_score:
            score = score.value
            return query, subject, score
            
        return query, subject

    return comprise 


@cAlgorithm
def NeedlemanWunsch_c(cFile_path):
    algorithm = cdll.LoadLibrary(cFile_path)
    return algorithm.NeedlemanWunsch

@cAlgorithm
def SmithWaterman_c(cFile_path):
    algorithm = cdll.LoadLibrary(cFile_path)
    return algorithm.SmithWaterman