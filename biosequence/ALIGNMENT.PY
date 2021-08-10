MATCH = 3
MISMATCH = -3
GAP = -3


class MatrixNode():
    def __init__(self,score):
        self.score = score
        self.up = False
        self.left = False
        self.upLeft = False

def SmithWaterman(seq1, seq2):
    matrix = [[MatrixNode(0) for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        matrix[i][0].up = True
    for j in range(1,len(seq2)+1):
        matrix[0][j].left = True
    maxNode = MatrixNode(0)
    # 计算得分矩阵
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            score = MATCH if seq1[i-1]==seq2[j-1] else MISMATCH
            up_score = matrix[i-1][j].score + GAP
            left_score = matrix[i][j-1].score + GAP
            upLeft_score = matrix[i-1][j-1].score + score
            matrix[i][j].score = max(up_score, left_score, upLeft_score,0)
            if matrix[i][j].score>=maxNode.score:
                maxNode = matrix[i][j]
                max_i = i
                max_j = j
            if matrix[i][j].score == up_score:
                matrix[i][j].up = True
            if matrix[i][j].score == left_score:
                matrix[i][j].left = True
            if matrix[i][j].score == upLeft_score:
                matrix[i][j].upLeft = True
    align_seq1 = ""
    align_seq2 = ""
    print(f"Max Score: {matrix[max_i][max_j].score}")
    while(matrix[max_i][max_j].score):
        node = matrix[max_i][max_j]
        if node.upLeft:
            align_seq1 += seq1[max_i-1]
            align_seq2 += seq2[max_j-1]
            max_i -= 1
            max_j -= 1
        elif node.up:
            align_seq1 += seq1[max_i-1]
            align_seq2 += "-"
            max_i -= 1
        elif node.left:
            align_seq1 += "-"
            align_seq2 += seq2[max_j-1]
            max_j -= 1
        else:
            break
    if max_i > max_j:
        align_seq2 += (max_i-max_j)*"-"
    if max_j > max_i:
        align_seq1 += (max_j-max_i)*"-"
    
    print(align_seq1[::-1])
    print(align_seq2[::-1])
    return

def NeedlemanWunsch(seq1, seq2):
    matrix = [[MatrixNode(0) for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        matrix[i][0].score = matrix[i-1][0].score + GAP
        matrix[i][0].up = True
    for j in range(1,len(seq2)+1):
        matrix[0][j].score = matrix[0][j-1].score + GAP
        matrix[0][j].left = True
    # 计算得分矩阵
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            score = MATCH if seq1[i-1]==seq2[j-1] else MISMATCH
            up_score = matrix[i-1][j].score + GAP
            left_score = matrix[i][j-1].score + GAP
            upLeft_score = matrix[i-1][j-1].score + score
            matrix[i][j].score = max(up_score, left_score, upLeft_score)
            if matrix[i][j].score == up_score:
                matrix[i][j].up = True
            if matrix[i][j].score == left_score:
                matrix[i][j].left = True
            if matrix[i][j].score == upLeft_score:
                matrix[i][j].upLeft = True
    align_seq1 = ""
    align_seq2 = ""
    while(i and j):
        node = matrix[i][j]
        if node.upLeft:
            align_seq1 += seq1[i-1]
            align_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif node.up:
            align_seq1 += seq1[i-1]
            align_seq2 += "-"
            i -= 1
        elif node.left:
            align_seq1 += "-"
            align_seq2 += seq2[j-1]
            j -= 1
    print(f"Max Score: {matrix[-1][-1].score}")
    print(align_seq1[::-1])
    print(align_seq2[::-1])



    

    return


if __name__ == '__main__':
    seq1 = "TGTTACGG"     # GTT_AC
    seq2 = "GGTTGACTA"    # GTTGAC
    SmithWaterman(seq1, seq2)

    seq1 = "GGATCGA"
    seq2 = "GAATTCAGTTA"
    NeedlemanWunsch(seq1, seq2)

