from config import AlignmentConfig

MATCH = AlignmentConfig.MATCH
MISMATCH = AlignmentConfig.MISMATCH
GAP_OPEN = AlignmentConfig.GAP_OPEN
GAP_EXTEND = AlignmentConfig.GAP_EXTEND


class MatrixNode:
    def __init__(self, score):
        self.score = score
        self.up = 0
        self.left = 0
        self.upLeft = 0

    def recordSource(self, up_score, left_score, upLeft_score):
        if self.score == up_score:
            self.up = 1
        if self.score == left_score:
            self.left = 1
        if self.score == upLeft_score:
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
    matrix = [[MatrixNode(0) for _ in range(length2 + 1)] for __ in range(length1 + 1)]

    for i in range(1, length1 + 1):
        matrix[i][0].score = matrix[i - 1][0].score + GAP
        matrix[i][0].up = True
    for j in range(1, length2 + 1):
        matrix[0][j].score = matrix[0][j - 1].score + GAP
        matrix[0][j].left = True

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

    up_score = matrix[current_i - 1][current_j].score + (
        GAP_EXTEND if matrix[current_i - 1][current_j].up else GAP_OPEN
    )
    left_score = matrix[current_i][current_j - 1].score + (
        GAP_EXTEND if matrix[current_i][current_j - 1].left else GAP_OPEN
    )
    upLeft_score = matrix[current_i - 1][current_j - 1].score + (
        MATCH if current_base1 == current_base2 else MISMATCH)

    return up_score, left_score, upLeft_score


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
    # 向左上移动说明seq1(column)的当前碱基与seq2(row)的当前碱基匹配
    if current_node.upLeft:
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
    # 向上移动说明seq2(row)的这个碱基要去匹配seq1(column)的下一个碱基
    # 即此时seq1的碱基匹配到的seq2为一个空位，因此seq2要在此增加一个‘-’
    # 即匹配的序列中，align_seq2[j]='-'，align_seq1[i]=seq1[i]
    elif current_node.up:
        sequence2 = "".join(
            [sequence2[: current_position[1]], "-", sequence2[current_position[1] :]]
        )
        current_position[0] -= 1

    return sequence1, sequence2


def SmithWaterman(query, subject):
    # initial score matrix
    matrix = initMatrix(len(query), len(subject), 0)
    maxNode = MatrixNode(0)

    # calculate score matrix
    for i in range(1, len(query) + 1):
        for j in range(1, len(subject) + 1):

            up_score, left_score, upLeft_score = getScore(
                i, j, matrix, query[i - 1], subject[j - 1]
            )

            matrix[i][j].score = max(up_score, left_score, upLeft_score, 0)
            matrix[i][j].recordSource(up_score, left_score, upLeft_score)

            if matrix[i][j].score >= maxNode.score:
                maxNode = matrix[i][j]
                position = [i, j]

    # Back Tracking
    while (node := matrix[position[0]][position[1]]).score:
        query, subject = backTracking(position, node, query, subject)

    # Align the sequences
    if position[0] > position[1]:
        subject = (position[0] - position[1]) * "." + subject
    if position[1] > position[0]:
        query = (position[1] - position[0]) * "." + query
    if len(query) < len(subject):
        query += "." * (len(subject) - len(query))
    if len(subject) < len(query):
        subject += "." * (len(query) - len(subject))

    return query, subject, maxNode.score, matrix


def NeedlemanWunsch(query, subject):
    # initial score matrix
    matrix = initMatrix(len(query), len(subject), 0)

    # calculate score matrix
    for i in range(1, len(query) + 1):
        for j in range(1, len(subject) + 1):

            up_score, left_score, upLeft_score = getScore(
                i, j, matrix, query[i - 1], subject[j - 1]
            )

            matrix[i][j].score = max(up_score, left_score, upLeft_score)
            matrix[i][j].recordSource(up_score, left_score, upLeft_score)

    # Back Tracking
    position = [i, j]
    while position[0] or position[1]:
        node = matrix[position[0]][position[1]]
        query, subject = backTracking(position, node, query, subject)

    return query, subject, matrix[-1][-1].score, matrix


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

    aligned_seq1, aligned_seq2, max_score, matrix = func(sequence1, sequence2)
    print(f"Max score: {max_score}")
    print(f"Aligned sequence1: {aligned_seq1}")
    print(f"Aligned sequence2: {aligned_seq2}")

    if show_matrix:
        print_matrix(matrix, sequence1, sequence1)
    print()


if __name__ == "__main__":
    seq1 = "TGTTACGG"  # TGTT_ACGG
    seq2 = "GGTTGACTA"  # GGTTGACTZ
    test(SmithWaterman, seq1, seq2)

    seq1 = "GGATCGA"  # GGA-TC-G--A
    seq2 = "GAATTCAGTTA"  # GAATTCAGTTA
    test(NeedlemanWunsch, seq1, seq2)

