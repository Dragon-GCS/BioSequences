#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GAP_CHAR '-';



typedef struct {
    int up;   
    int left;   
    int upLeft;   
    int score;
}matrixNode, *mNode;

mNode** initMatrix(int row_length, int column_length, int gap_open, int gap_extend);

int max3(int a, int b, int c);

void backTracking(mNode node, int* current_i, int* current_j, char query_base, char subject_base, char* align_query, char* align_subject, int index);
void reverseStr(char* str);
void release(mNode** matrix, int rows, int columns);

void assignAlignParameter(char gap_char, int match, int mismatch, int gap_open, int gap_extend);
void NeedlemanWunsch(char* query, char* subject, 
                     char* aligned_query, char* aligned_subject,
                     int* score,  
                     int match, int mismatch, 
                     int gap_open, int gap_extend);

void SmithWaterman(char* query, char* subject, 
                   char* aligned_query, char* aligned_subject, 
                   int* score,  
                   int match, int mismatch, 
                   int gap_open, int gap_extend);


mNode** initMatrix(int row_length, int column_length, int gap_open, int gap_extend)
{   
    int i=0, j=0;
    mNode **matrix;

    if ((matrix = (mNode **) malloc(sizeof(mNode*) * (row_length + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
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
            matrix[i][j]->up = 0;
            matrix[i][j]->left = 0;
            matrix[i][j]->upLeft = 0;
        }
    }
    matrix[0][0]->score = 0;
    for (i = 1; i < row_length; i++) {
        matrix[i][0]->score = gap_open + gap_extend * (i - 1);
        matrix[i][0]->up = 1;
    }
    for (j = 1; j < column_length; j++) {
        matrix[0][j]->score = gap_open + gap_extend * (i - 1);
        matrix[0][j]->left = 1;
    }
    return matrix;
}


int max3(int a, int b, int c) 
{
    int f = a > b ? a : b;
    return f > c ? f : c;
}


void backTracking(mNode node, int* current_i, int* current_j, char query_base, char subject_base, char* align_query, char* align_subject, int index)
{
    if( node->upLeft )
    {
        *current_i -= 1;
        *current_j -= 1;
        align_query[index] = query_base;
        align_subject[index] = subject_base;
    }
    else if( node->left )
    {
        *current_j -= 1;
        align_query[index] = GAP_CHAR;
        align_subject[index] = subject_base;
    }
    else if( node->up )
    {
        *current_i -= 1;
        align_query[index] = query_base;
        align_subject[index] = GAP_CHAR;
    }
}


void release(mNode** matrix, int rows, int columns)
{   
    for (int i = 0; i <= rows; i++) {
        for (int j = 0; j <= columns; j++) {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}


void reverseStr(char* str)
{
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
                     int* score,  
                     int match, int mismatch, 
                     int gap_open, int gap_extend) 
{
    int i=0, j=0, index = 0;
    int up_score, upLeft_score, left_score, node_score;

    int rows = strlen(query);
    int columns = strlen(subject);

    mNode **score_matrix = initMatrix(rows+1, columns+1, gap_open, gap_extend);
    
    // Calculate the score
    for (i = 1; i <= rows; i++) {
        for (j = 1; j <= columns; j++) {

            up_score = score_matrix[i - 1][j]->score + ( score_matrix[i - 1][j]->up ? gap_extend : gap_open);
            left_score = score_matrix[i][j -1]->score + (score_matrix[i][j -1]->left ? gap_extend : gap_open);
            upLeft_score = score_matrix[i - 1][j -1]->score + (query[i-1] == subject[j-1] ? match : mismatch);

            node_score = max3(up_score, left_score, upLeft_score);

            score_matrix[i][j]->score = node_score;
            if (up_score == node_score) score_matrix[i][j]->up = 1;
            if (left_score == node_score) score_matrix[i][j]->left = 1;
            if (upLeft_score == node_score) score_matrix[i][j]->upLeft = 1;
        }
    }

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
                   int* score,  
                   int match, int mismatch, 
                   int gap_open, int gap_extend)
{
    int i=0, j=0, index = 0;
    int up_score, upLeft_score, left_score, node_score;
    int max_i, max_j;
    *score = 0;

    int rows = strlen(query);
    int columns = strlen(subject);

    mNode **score_matrix = initMatrix(rows+1, columns+1, 0, 0);
    
    // Calculate the score
    for (i = 1; i <= rows; i++) {
        for (j = 1; j <= columns; j++) {

            up_score = score_matrix[i - 1][j]->score + ( score_matrix[i - 1][j]->up ? gap_extend : gap_open);
            left_score = score_matrix[i][j -1]->score + (score_matrix[i][j -1]->left ? gap_extend : gap_open);
            upLeft_score = score_matrix[i - 1][j -1]->score + (query[i-1] == subject[j-1] ? match : mismatch);

            node_score = max3(up_score, left_score, upLeft_score);
            if (node_score < 0)  node_score = 0;

            score_matrix[i][j]->score = node_score;
            if (up_score == node_score) score_matrix[i][j]->up = 1;
            if (left_score == node_score) score_matrix[i][j]->left = 1;
            if (upLeft_score == node_score) score_matrix[i][j]->upLeft = 1;

            if (node_score >= (*score))
            {   
                *score = node_score;
                max_i = i; max_j = j;
            }
            
        }
    }

    
    // get the aligned sequence by back-tracking
    i = max_i; j=max_j;
    while(score_matrix[i][j]->score)
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
