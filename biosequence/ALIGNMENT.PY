MATCH = 10
MISMATCH = -2
GAP = -50

class MatrixNode():
    def __init__(self,score):
        self.score = score
        self.up = False
        self.left = False
        self.upLeft = False


def SmithWaterman(seq1, seq2):
    # initial score matrix
    matrix = [[MatrixNode(0) for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        matrix[i][0].up = True
    for j in range(1,len(seq2)+1):
        matrix[0][j].left = True
    maxNode = MatrixNode(0)
    # calculate score matrix
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

    print(f"Max Score: {matrix[max_i][max_j].score}")
    # backward
    while(matrix[max_i][max_j].score):
        node = matrix[max_i][max_j]
        if node.upLeft:
            max_i -= 1
            max_j -= 1
        elif node.up:
            seq2 = "".join([seq2[:max_j],"-",seq2[max_j:]]) 
            max_i -= 1
        elif node.left:
            seq1 = "".join([seq1[:max_i],"-",seq1[max_i:]]) 
            max_j -= 1
        
    if max_i > max_j:
        seq2 = (max_i-max_j)*"-" + seq2
    if max_j > max_i:
        seq1 = (max_j-max_i)*"-" + seq1
    if len(seq1) < len(seq2):
        seq1 += "-"*(len(seq2)-len(seq1))
    if len(seq2) < len(seq1):
        seq2 += "-"*(len(seq1)-len(seq2))
    
    print(seq1)
    print(seq2)
    return matrix


def NeedlemanWunsch(seq1, seq2):
    # initial score matrix
    matrix = [[MatrixNode(0) for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        matrix[i][0].score = matrix[i-1][0].score + GAP
        matrix[i][0].up = True
    for j in range(1,len(seq2)+1):
        matrix[0][j].score = matrix[0][j-1].score + GAP
        matrix[0][j].left = True
    # calculate score matrix
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

    while(i or j):
        node = matrix[i][j]
        # 向左上移动说明seq1(column)的当前碱基与seq2(row)的当前碱基匹配
        if node.upLeft:
            i -= 1
            j -= 1
        # 向上移动说明seq2(row)的这个碱基要去匹配seq1(column)的下一个碱基
        # 即此时seq1的碱基匹配到的seq2为一个空位，因此seq2要在此增加一个‘-’
        # 即匹配的序列中，align_seq2[j]='-'，align_seq1[i]=seq1[i]
        elif node.up:
            seq2 = "".join([seq2[:j],"-",seq2[j:]]) 
            i -= 1
        # 向左移动说明seq1(column)的这个碱基要去匹配seq2(row)的下一个碱基
        # 即此时seq2的碱基匹配到的seq1碱基为一个空位，因此seq1要在此增加一个‘-’
        # 即匹配的序列中，align_seq1[i]='-'，align_seq2[j]=seq2[j]
        elif node.left:
            seq1 = "".join([seq1[:i],"-",seq1[i:]]) 
            j -= 1
    print(f"Max Score: {matrix[-1][-1].score}")
    print(seq1)
    print(seq2)
    return matrix


def print_matrix(matrix, seq1, seq2):
    print(''.join(f"{bp:4s}".center(4) for bp in ("  "+ seq2)))
    seq1 = " " + seq1
    for i in range(len(seq1)):
        print(seq1[i],end="  ")
        print(" ".join(list(f"{str(node.score):3s}".center(3) for node in matrix[i])))


def test(func):
    from random import choice, randint

    def generate(num = 0):
        if not num:
            num = randint(1,100)
        return ''.join([choice(["A","T", "C", "G"]) for i in range(num)])

    seq1 = generate(5)
    seq2 = generate(10)
    print(f"\nSequence 1:{seq1}")
    print(f"Sequence 2:{seq2}")
    matrix = func(seq1, seq2)
    print_matrix(matrix, seq1, seq2)
        

if __name__ == '__main__':
    seq1 = "TGTTACGG"           # TGTT_ACGG
    seq2 = "GGTTGACTA"          # GGTTGACTZ
    SmithWaterman(seq1, seq2)

    seq1 = "GGATCGA"            # GGA-TC-G--A
    seq2 = "GAATTCAGTTA"        # GAATTCAGTTA
    NeedlemanWunsch(seq1, seq2)
    
    test(NeedlemanWunsch)
    test(SmithWaterman)


