#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algorithm.h"
#define GAP_CHAR '-';

mNode** initMatrix(int row_length, int column_length)
{   /*
     *Description:  Initial the score matrix
     *Input:    
        @row_length:    The num of rows, equals to the len(query) + 1
        @column_length: The num of columns, equals to the len(subject) + 1
     *Output: None
     *Return:
        @matrix         The pointer of score matrix
     */
    int i=0, j=0;
    mNode **matrix;

    if ((matrix = (mNode **) malloc(sizeof(mNode*) * (row_length + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    /* create row_length x columns_length nodes*/ 
    for (i = 0; i < row_length; i++) {
        if ((matrix[i] = (mNode *) malloc(sizeof(mNode) * (column_length + 1))) == NULL)
        {
            fputs("Error: Out of space!\n", stderr);
            exit(1);     
        }
        for (j = 0; j < column_length; j++) 
        {
            if ((matrix[i][j] = (mNode) malloc(sizeof(matrixNode))) == NULL) 
            {
                fputs("Error: Out of space!\n", stderr);
                exit(1);     
            }
            matrix[i][j]->uScore = 0; matrix[i][j]->lScore = 0; matrix[i][j]->mScore = 0;
            matrix[i][j]->up = 0;     matrix[i][j]->left = 0;   matrix[i][j]->upLeft = 0;
            matrix[i][j]->score = 0;
        }
    }
    return matrix;
}

float max2(float a, float b) 
{   /*
     *Description: return the maxium of two value
     */
    return a > b ? a : b;
}

float max3(float a, float b, float c) 
{   /*
     *Description: return the maxium of three value
     */
    float f = a > b ? a : b;
    return f > c ? f : c;
}


void backTracking(mNode node, int* current_i, int* current_j, char query_base, char subject_base, char* align_query, char* align_subject, int index)
{   /*
     *Description: Process the back-tracking
     *Input:
        @node:          current node
        @current_i:     if node->up = 1, current_i--
        @current_j:     if node->left = 1, current_j--
        @query_base:    base in query to compare
        @subject_base:  base in subject to compare
        @index:         index for aligned string
     *Output:
        @align_query:   aligned sting by reverse order
        @align_subject: aligned sting by reverse order
     *Return:           None
     */
    if( node->left )
    {
        (*current_j)--;
        align_query[index] = GAP_CHAR;
        align_subject[index] = subject_base;
    }
   else if( node->upLeft )
    {
        (*current_i)--;
        (*current_j)--;
        align_query[index] = query_base;
        align_subject[index] = subject_base;
    } 
    else if( node->up )
    {
        (*current_i)--;
        align_query[index] = query_base;
        align_subject[index] = GAP_CHAR;
    } 
}


void release(mNode** matrix, int rows, int columns)
{   /*
     *Description: free the matrix
     *Input:
        @matrix:    the pointer of matrix
        @rows:      Rows num of matrix
        @columns    Columns of matrix
     */
    for (int i = 0; i <= rows; i++) {
        for (int j = 0; j <= columns; j++) {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}


void reverseStr(char* str)
{   /*
     *Description:  Reverse a string
     *Input:     
        @str:   Pointer of a string
     */
    int length = strlen(str);
    int pair = length / 2;
    int last;
    char temp;
    for( int i=0; i<pair; i++)
    {
        last = length - i - 1;
        temp = str[i]; str[i] = str[last]; str[last] = temp;
    }
}


void NeedlemanWunsch(char* query, char* subject, 
                     char* aligned_query, char* aligned_subject, 
                     float* score,  
                     float match, float mismatch, 
                     float gap_open, float gap_extend) 
{   /*
     *Description:  Align query and subject by Needleman-Wunsch
     *Input:
        @query:     Sequence string 1
        @subject:   Sequence string 2
        @match:     Score when two base are same
        @mismatch:  Score when two base are different
        @gap_open:  Score when gap appear
        @gap_extend:Score when gap extend
     *Output:
        @aligned_query:     Sequence 1 after aligned
        @aligned_subject:   Sequence 2 after aligned
        @score:             Max score of align
     *Return:   None
     */
    float match_score;
    int i = 0, j = 0, index = 0;
    int rows = strlen(query);
    int columns = strlen(subject);
    
    // Initial the score matrix
    mNode **score_matrix = initMatrix(rows+1, columns+1);

    // Assign init score to the 1st column's and 1st row's node
    score_matrix[0][0]->uScore = gap_open;
    score_matrix[0][0]->lScore = gap_open;
    // Node's score = maxium of uScore, lScore and mScore
    // Beacause gap open is minus and mScore is zero, the node's score shoud be zero

    for (i = 1; i <= rows; i++) {
        score_matrix[i][0]->uScore = 2 * gap_open + gap_extend * (i - 1);
        score_matrix[i][0]->lScore = gap_open + gap_extend * (i - 1);
        score_matrix[i][0]->mScore = gap_open + gap_extend * (i - 1);
        score_matrix[i][0]->score = max3(score_matrix[i][0]->uScore, 
                                         score_matrix[i][0]->lScore, 
                                         score_matrix[i][0]->mScore);
        score_matrix[i][0]->up = 1;
    }
    for (j = 1; j <= columns; j++) {
        score_matrix[0][j]->uScore = gap_open + gap_extend * (i - 1);
        score_matrix[0][j]->lScore = 2 * gap_open + gap_extend * (i - 1);
        score_matrix[0][j]->mScore = gap_open + gap_extend * (i - 1);
        score_matrix[0][j]->score = max3(score_matrix[0][j]->uScore, 
                                         score_matrix[0][j]->lScore, 
                                         score_matrix[0][j]->mScore);
        score_matrix[0][j]->left = 1;
    }

    // Calculate the score
    for (i = 1; i <= rows; i++) 
    {
        for (j = 1; j <= columns; j++) 
        {
            score_matrix[i][j]->uScore = max2(score_matrix[i - 1][j]->mScore + gap_open,
                                              score_matrix[i - 1][j]->uScore + gap_extend);
            score_matrix[i][j]->lScore = max2(score_matrix[i][j - 1]->mScore + gap_open,
                                              score_matrix[i][j - 1]->lScore + gap_extend);

            match_score = query[i-1] == subject[j-1] ? match : mismatch;
            score_matrix[i][j]->mScore = score_matrix[i - 1][j - 1]->score + match_score;
            

            score_matrix[i][j]->score = max3(score_matrix[i][j]->uScore, 
                                             score_matrix[i][j]->lScore, 
                                             score_matrix[i][j]->mScore);

            if (score_matrix[i][j]->uScore == score_matrix[i][j]->score) 
            {
                score_matrix[i][j]->up = 1;
                if (i > 1)
                {
                if ((score_matrix[i][j]->score == score_matrix[i - 1][j]->uScore + gap_extend) & (score_matrix[i - 1][j] ->score == score_matrix[i - 1][j]->uScore))
                    {score_matrix[i - 1][j]->up=0; score_matrix[i - 1][j]->upLeft = 0; score_matrix[i - 1][j]->up = 1;}
                if ((score_matrix[i][j]->score == score_matrix[i - 1][j]->mScore + gap_open) & (score_matrix[i - 1][j] ->score == score_matrix[i - 1][j]->mScore))
                    {score_matrix[i - 1][j]->up=0; score_matrix[i - 1][j]->upLeft = 1; score_matrix[i - 1][j]->up = 0;}
                } 
            }
            if (score_matrix[i][j]->lScore == score_matrix[i][j]->score) 
            {
                score_matrix[i][j]->left = 1;
                if (i > 1)
                {
                if ((score_matrix[i][j]->score == score_matrix[i][j - 1]->lScore + gap_extend) & (score_matrix[i][j - 1] ->score == score_matrix[i][j - 1]->lScore))
                    {score_matrix[i][j-1]->left = 1; score_matrix[i][j-1]->upLeft = 0; score_matrix[i][j-1]->up = 0;}
                if ((score_matrix[i][j]->score == score_matrix[i][j - 1]->mScore + gap_open) & (score_matrix[i][j - 1] ->score == score_matrix[i][j - 1]->mScore))
                    {score_matrix[i][j-1]->left = 0; score_matrix[i][j-1]->upLeft = 1; score_matrix[i][j-1]->up = 0;}
                }
                
            }
            if (score_matrix[i][j]->mScore == score_matrix[i][j]->score) score_matrix[i][j]->upLeft = 1;
        }
    }
    // record the max score
    *score = score_matrix[rows][columns]->score;
    
    // get the aligned sequence by back-tracking
    i = rows; j=columns;
    while((i || j) )
    {
        backTracking(score_matrix[i][j], &i, &j, query[i-1], subject[j-1], aligned_query, aligned_subject, index);
        index++;
    }
    aligned_query[index] = '\0';
    aligned_subject[index] = '\0';

    reverseStr(aligned_query);
    reverseStr(aligned_subject);
    
    release(score_matrix, rows, columns);
}


void SmithWaterman(char* query, char* subject, 
                   char* aligned_query, char* aligned_subject, 
                   float* score,  
                   float match, float mismatch, 
                   float gap_open, float gap_extend)
{   /*
     *Description:  Align query and subject by Smith-Waterman
     *Input:
        @query:     Sequence string 1
        @subject:   Sequence string 2
        @match:     Score when two base are same
        @mismatch:  Score when two base are different
        @gap_open:  Score when gap appear
        @gap_extend:Score when gap extend
     *Output:
        @aligned_query:     Sequence 1 after aligned
        @aligned_subject:   Sequence 2 after aligned
        @score:             Max score of align
     *Return:   None
     */
    float match_score;
    int i = 0, j = 0, index = 0;
    int max_i = 0, max_j = 0;
    *score = 0;

    int rows = strlen(query);
    int columns = strlen(subject);

    // Initial the score matrix
    mNode **score_matrix = initMatrix(rows+1, columns+1);

    // Calculate the score
    for (i = 1; i <= rows; i++) 
    {
        for (j = 1; j <= columns; j++)
        {

            score_matrix[i][j]->uScore = max3(score_matrix[i - 1][j]->mScore + gap_open,
                                              score_matrix[i - 1][j]->uScore + gap_extend,
                                              0);
            score_matrix[i][j]->lScore = max3(score_matrix[i][j - 1]->mScore + gap_open,
                                              score_matrix[i][j - 1]->lScore + gap_extend,
                                              0);

            match_score = query[i-1] == subject[j-1] ? match : mismatch;
            score_matrix[i][j]->mScore = max2(0, score_matrix[i - 1][j - 1]->score + match_score);
            

            score_matrix[i][j]->score = max3(score_matrix[i][j]->uScore, 
                                             score_matrix[i][j]->lScore, 
                                             score_matrix[i][j]->mScore);

            if( score_matrix[i][j]->score != 0)
            {      
                if (score_matrix[i][j]->uScore == score_matrix[i][j]->score) 
                {
                    score_matrix[i][j]->up = 1;
                    if (i > 1)
                    {
                    if ((score_matrix[i][j]->score == score_matrix[i - 1][j]->uScore + gap_extend) & (score_matrix[i - 1][j] ->score == score_matrix[i - 1][j]->uScore))
                        {score_matrix[i - 1][j]->up=0; score_matrix[i - 1][j]->upLeft = 0; score_matrix[i - 1][j]->up = 1;}
                    if ((score_matrix[i][j]->score == score_matrix[i - 1][j]->mScore + gap_open) & (score_matrix[i - 1][j] ->score == score_matrix[i - 1][j]->mScore))
                        {score_matrix[i - 1][j]->up=0; score_matrix[i - 1][j]->upLeft = 1; score_matrix[i - 1][j]->up = 0;}
                    }
                }
                if (score_matrix[i][j]->lScore == score_matrix[i][j]->score) 
                {
                    score_matrix[i][j]->left = 1;
                    if (i > 1)
                    {
                    if ((score_matrix[i][j]->score == score_matrix[i][j - 1]->lScore + gap_extend) & (score_matrix[i][j - 1] ->score == score_matrix[i][j - 1]->lScore))
                        {score_matrix[i][j-1]->left = 1; score_matrix[i][j-1]->upLeft = 0; score_matrix[i][j-1]->up = 0;}
                    if ((score_matrix[i][j]->score == score_matrix[i][j - 1]->mScore + gap_open) & (score_matrix[i][j - 1] ->score == score_matrix[i][j - 1]->mScore))
                        {score_matrix[i][j-1]->left = 0; score_matrix[i][j-1]->upLeft = 1; score_matrix[i][j-1]->up = 0;}
                    }
                    
                }      
                if (score_matrix[i][j]->mScore == score_matrix[i][j]->score) score_matrix[i][j]->upLeft = 1;
                
                // Recode the max score and its positon
                if (score_matrix[i][j]->score >= (*score) )
                {   
                    *score = score_matrix[i][j]->score;
                    max_i = i; max_j = j;
                }
            }
        }
    }

    
    // get the aligned sequence by back-tracking
    i = max_i; j=max_j;
    while(score_matrix[i][j]->score != 0)
    {
        backTracking(score_matrix[i][j], &i, &j, query[i-1], subject[j-1], aligned_query, aligned_subject, index);
        index++;
    }
    aligned_query[index] = '\0';
    aligned_subject[index] = '\0';

    reverseStr(aligned_query);
    reverseStr(aligned_subject);
    
    release(score_matrix, rows, columns);
}