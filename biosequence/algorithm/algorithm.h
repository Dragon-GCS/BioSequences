
/* Unit in Score Matrix */
typedef struct {
    int up;         /* whether go up node when back-tracking */
    int left;       /* whether go left node when back-tracking */
    int upLeft;     /* whether go up-left node when back-tracking */
    float score;    /* node's final score */
    float uScore;   /* score if route from up node */
    float lScore;   /* score if route from left node */
    float mScore;   /* score if route from up-left node */
}matrixNode, *mNode;

mNode** initMatrix(int row_length, int column_length);
float max3(float a, float b, float c);
void backTracking(mNode node, int* current_i, int* current_j, char query_base, char subject_base, char* align_query, char* align_subject, int index);
void reverseStr(char* str);
void release(mNode** matrix, int rows, int columns);

void NeedlemanWunsch(char* query, char* subject, 
                     char* aligned_query, char* aligned_subject,
                     float* score,  
                     float match, float mismatch, 
                     float gap_open, float gap_extend);

void SmithWaterman(char* query, char* subject, 
                   char* aligned_query, char* aligned_subject, 
                   float* score,  
                   float match, float mismatch, 
                   float gap_open, float gap_extend);

