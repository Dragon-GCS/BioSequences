from biosequence.alignment._align_utils import *


def SmithWaterman(query, subject, return_score=False):
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
    if return_score:
        return query, subject, maxNode.score
            
    return query, subject
    


def NeedlemanWunsch(query, subject, return_score=False):
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
    if return_score:
        return query, subject, matrix[-1][-1].score
            
    return query, subject
    

@cAlgorithm
def NeedlemanWunsch_c(cFile_path):
    algorithm = cdll.LoadLibrary(cFile_path)
    return algorithm.NeedlemanWunsch

@cAlgorithm
def SmithWaterman_c(cFile_path):
    algorithm = cdll.LoadLibrary(cFile_path)
    return algorithm.SmithWaterman