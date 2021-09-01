import platform     # detect system for using .dll or .so to boost alignment
import os

from functools import wraps
from ctypes import *
from biosequence.config import AlignmentConfig


class MatrixNode:
    def __init__(self, lScore, uScore, mScore):
        
        self.uScore = 0
        self.lScore = 0
        self.mScore = 0
        self.score = max(lScore,uScore,mScore)
        self.up = 0
        self.left = 0
        self.upLeft = 0

    def getScore(self):
        self.score = max(self.lScore,self.uScore,self.mScore)

    def recordSource(self):
        if self.score == self.uScore:
            self.up = 1
        if self.score == self.lScore:
            self.left = 1
        if self.score == self.mScore:
            self.upLeft = 1


def initMatrix(length1, length2, GAP):
    """
    Initialization the score matrix
    Args:
        length1: The length of sequence1 as the num of rows
        length2: The length of sequence2 as the num of columns
        GAP: GAP score for SmithWaterman
    Returns:
        matrix: Initiallized matrix
    """
    matrix = [[MatrixNode(0,0,0) for _ in range(length2 + 1)] for __ in range(length1 + 1)]
    matrix[0][0].lScore = AlignmentConfig.GAP_OPEN
    matrix[0][0].uScore = AlignmentConfig.GAP_OPEN
    matrix[0][0].mScore = 0
    matrix[0][0].getScore()
    
    for i in range(1, length1 + 1):
        matrix[i][0].uScore = 2 * AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[i][0].lScore = AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[i][0].mScore = AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[i][0].getScore()
        matrix[i][0].up = 1

    for j in range(1, length2 + 1):
        matrix[0][j].uScore = AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[0][j].lScore = 2 * AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[0][j].mScore = AlignmentConfig.GAP_OPEN + (i - 1) * AlignmentConfig.GAP_EXTEND
        matrix[0][j].getScore()
        matrix[0][j].left = 1
        
    return matrix


def getScore(current_i, current_j, matrix, current_base1, current_base2):
    """
    Calculate the score of current node
    Args:
        current_i: Index of current row
        current_j: Index of current column
        matrix: The score matrix
        current_base1: The pairing base of sequence1
        current_base2: The pairing base of sequence2
    Returns:
        up_score: The score if the route came from up node
        left_score: The score if the route came from left node
        upleft_score: The score if the route came from upleft node
    """
    match_score = AlignmentConfig.MATCH if current_base1 == current_base2 else AlignmentConfig.MISMATCH
    
    uScore = max(
                matrix[current_i - 1][current_j].mScore + AlignmentConfig.GAP_OPEN, 
                matrix[current_i - 1][current_j].uScore + AlignmentConfig.GAP_EXTEND)

    lScore = max(
                matrix[current_i][current_j - 1].mScore + AlignmentConfig.GAP_OPEN,
                matrix[current_i][current_j - 1].mScore + AlignmentConfig.GAP_EXTEND)

    mScore = max(
                matrix[current_i - 1][current_j - 1].mScore + match_score,
                matrix[current_i - 1][current_j - 1].uScore + match_score,
                matrix[current_i - 1][current_j - 1].lScore + match_score)
    return lScore, uScore, mScore
    


def backTracking(current_position, current_node, sequence1, sequence2):
    """
    Back tracking from the current position
    Args:
        current_position: the position list which is tracking now, [row_index, column_index]
        current_node: the node of score matrix which is tracking now
        sequence1: sequence1
        sequence2: sequence2 
    Returns:
        aligned_sequence1: aligned sequence1
        aligned_sequence2: aligned sequence2       
    """
    # 向上移动说明seq2(row)的这个碱基要去匹配seq1(column)的下一个碱基
    # 即此时seq1的碱基匹配到的seq2为一个空位，因此seq2要在此增加一个‘-’
    # 即匹配的序列中，align_seq2[j]='-'，align_seq1[i]=seq1[i]
    if current_node.up:
        sequence2 = "".join(
            [sequence2[: current_position[1]], "-", sequence2[current_position[1] :]]
        )
        current_position[0] -= 1

    # 向左上移动说明seq1(column)的当前碱基与seq2(row)的当前碱基匹配
    elif current_node.upLeft:
        current_position[0] -= 1
        current_position[1] -= 1
    # 向左移动说明seq1(column)的这个碱基要去匹配seq2(row)的下一个碱基
    # 即此时seq2的碱基匹配到的seq1碱基为一个空位，因此seq1要在此增加一个‘-’
    # 即匹配的序列中，align_seq1[i]='-'，align_seq2[j]=seq2[j]
    elif current_node.left:
        sequence1 = "".join(
            [sequence1[: current_position[0]], "-", sequence1[current_position[0] :]]
        )
        current_position[1] -= 1
    return sequence1, sequence2

    


def test(func, sequence1="", sequence2="", show_matrix=False):
    from random import choice, randint

    def generate(num=0):
        if not num:
            num = randint(1, 100)
        return "".join([choice(["A", "T", "C", "G"]) for _ in range(num)])
    
    def print_matrix(matrix, seq1, seq2):
        print("".join(f"{bp:4s}".center(4) for bp in ("  " + seq2)))
        seq1 = " " + seq1
        for i in range(len(seq1)):
            print(seq1[i], end="  ")
            print(" ".join(list(f"{node.score:2d}".center(3) for node in matrix[i])))

    if not sequence1:
        sequence1 = generate(20)
    if not sequence2:
        sequence2 = generate(20)

    print(func.__name__ + ":")
    print(f"Sequence 1:{sequence1}")
    print(f"Sequence 2:{sequence2}")

    aligned_seq1, aligned_seq2, max_score = func(sequence1, sequence2)
    print(f"Max score: {max_score}")
    print(f"Aligned sequence1: {aligned_seq1}")
    print(f"Aligned sequence2: {aligned_seq2}")

    if show_matrix:
        print_matrix(matrix, sequence1, sequence1)
    print()


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