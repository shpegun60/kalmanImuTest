#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "smart_assert/smart_assert.h"

int showmat(Mat* A, char * name)
{
    if(A == NULL) {
        printf("matrix is not exist!!! \n\n");
        return MAT_FAIL;
    }

    printf("%s [%d:%d]\t\n", name, A->row, A->col);
    if((A->row == 0) && (A->col == 0)) {
        printf("[]");
        return MAT_SUCC;
    }

    char txt[100];
    printf("[");
    for(unsigned int i = 0; i < A->row; i++){
        for (unsigned int j = 0; j < A->col; j++) {
            sprintf(txt, "\t%.5f", (float)A->data[i][j]);
            printf("%s", txt);
        }
        if(i == (A->row - 1)) {
            printf("\t]\n");
        } else {
            printf("\n");
        }
    }
    printf("\n");
    return MAT_SUCC;
}


/* make a zero matrix of given dimensions */
Mat *matrixCreate(unsigned int r, unsigned int c)
{
    M_Assert_BreakSaveCheck((r == 0 || c == 0), "matrixCreate: Give me positive values for dimensions genius", return NULL);

    Mat *m;
    m = (Mat *)malloc(sizeof(Mat));
    M_Assert_BreakSaveCheck((m == NULL), "matrixCreate: no memories for allocation matrix", return NULL);

    m->row = r;
    m->col = c;

    m->data = (MAT_TYPE**)malloc(r * sizeof(MAT_TYPE*));
    M_Assert_BreakSaveCheck((m->data == NULL), "matrixCreate: no memories for allocation data", return NULL);
    for (unsigned int i = 0; i < r; ++i) {
        m->data[i] = (MAT_TYPE*)calloc(sizeof(MAT_TYPE), c);
        M_Assert_Break((m->data[i] == NULL), "matrixCreate: no memories for allocation", return NULL);
    }
    return m;
}

int matrixInitFromArr(Mat* A, MAT_TYPE* arr)
{
    M_Assert_Break((A == NULL || arr == NULL), "matrixInitFromArr: incorrect input values", return MAT_FAIL);
    for(unsigned int i = 0; i < A->row; ++i) {
        for(unsigned int j = 0; j < A->col; ++j) {
            A->data[i][j] =  *(arr + (i * A->col) + j);
        }
    }
    return MAT_SUCC;
}

int matrixInitFromArr_T(Mat* A, MAT_TYPE* arr)
{
    M_Assert_Break((A == NULL || arr == NULL), "matrixInitFromArr_T: incorrect input values", return MAT_FAIL);
    for(unsigned int i = 0; i < A->row; ++i) {
        for(unsigned int j = 0; j < A->col; ++j) {
            A->data[i][j] =  *(arr + (j * A->row) + i);
        }
    }
    return MAT_SUCC;
}



int getResultRawSize(Mat* A, Mat* B, unsigned int * resRaw, unsigned int * resCall)
{
    M_Assert_Break((!A || !B || !resRaw || !resCall), "getResultRawSize: incorrect input values", return MAT_FAIL);
    if(A->col != B->row) {
        *resRaw = 0;
        *resCall = 0;
        M_Assert_Break((A->col != B->row), "getResultRawSize: incorrect length`s", return MAT_FAIL);
        return MAT_FAIL;
    }
    *resRaw = A->row;
    *resCall = B->col;
    return MAT_SUCC;
}


Mat *createResultMulMatrix(Mat* A, Mat* B)
{
    M_Assert_Break((!A || !B), "createResultMulMatrix: incorrect input values", return NULL);
    unsigned int raw;
    unsigned int coll;

    if(getResultRawSize(A, B, &raw, &coll) != MAT_SUCC) {
        return NULL;
    }

    return matrixCreate(raw, coll);
}

Mat *createResultTransMatrix(Mat* A)
{
    M_Assert_Break((!A), "createResultTransMatrix: incorrect input values", return NULL);
    return matrixCreate(A->col, A->row);
}


/* free memory associated with the matrix  */
int destroy_matrix(Mat **m)
{
    M_Assert_BreakSaveCheck(((*m) == NULL || m == NULL), "destroy_matrix: incorrect input values", return MAT_FAIL);

    /* Code for further processing and free the
       dynamically allocated memory */

    for(unsigned int i = 0; i < (*m)->row; ++i) {
        free((*m)->data[i]);
    }

    free((*m)->data);
    free(*m);
    *m = NULL;
    return MAT_SUCC;
}



/* enter 1s along the main diagonal */
Mat *eye(unsigned int length)
{
    M_Assert_BreakSaveCheck((length == 0), "eye: Give me positive values for length genius", return NULL);

    Mat *m;
    m = matrixCreate(length, length);
    for(unsigned int i = 0; i < length; ++i){
        m->data[i][i] = (MAT_TYPE)1;
    }
    return m;
}

Mat* ones(unsigned int r, unsigned int c, MAT_TYPE d)
{
    M_Assert_BreakSaveCheck((r == 0 || c == 0), "ones: Give me positive values for dimensions genius", return NULL);

    Mat* M = matrixCreate(r,c);
    for(unsigned int i = 0; i < M->row; ++i){
        for(unsigned int j = 0; j < M->col; ++j){
            M->data[i][j] = d;
        }
    }
    return M;
}

int add(Mat* A, Mat* B, Mat* Dest) // hardness function: r1 * c1
{
    M_Assert_Break((!A || !B || !Dest), "add: incorrect input values", return MAT_FAIL);
    M_Assert_Break(( (A->row != B->row) || (A->col != B->col) || (A->row > Dest->row)  || (A->col > Dest->col) ), "add: incorrect length", return MAT_FAIL);

    unsigned int r = A->row;
    unsigned int c = A->col;

    for(unsigned int i = 0; i < r; ++i){
        for(unsigned int j = 0; j < c; ++j){
            Dest->data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
    return MAT_SUCC;
}

int sub(Mat* A, Mat* B, Mat* Dest) // hardness function: r1 * c1
{
    M_Assert_Break((!A || !B || !Dest), "sub: incorrect input values", return MAT_FAIL);
    M_Assert_Break(( (A->row != B->row) || (A->col != B->col) || (A->row > Dest->row)  || (A->col > Dest->col) ), "sub: incorrect length", return MAT_FAIL);
    unsigned int r = A->row;
    unsigned int c = A->col;

    for(unsigned int i = 0; i < r; ++i) {
        for(unsigned int j = 0; j < c; ++j){
            Dest->data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return MAT_SUCC;
}

int scalarmultiply(Mat* A, Mat* Dest, MAT_TYPE scalar)
{
    M_Assert_Break((!A || !Dest), "scalarmultiply: incorrect input values", return MAT_FAIL);
    M_Assert_Break((  (A->row > Dest->row)  || (A->col > Dest->col) ), "scalarmultiply: incorrect destinalion length", return MAT_FAIL);
    unsigned int r = A->row;
    unsigned int c = A->col;

    for(unsigned int i = 0; i < r; ++i) {
        for(unsigned int j = 0; j < c; ++j){
            Dest->data[i][j] = A->data[i][j] * scalar;
        }
    }
    return MAT_SUCC;
}

int multiply(Mat* A, Mat* B, Mat* Dest) // hardness function: 2 * r1 * c2 * r2
{    
    M_Assert_Break((!A || !B || !Dest), "multiply: incorrect input values", return MAT_FAIL);
    M_Assert_Break(  (          (A->col != B->row) && !((A->row == 1) && (A->col == 1)) && !((B->row == 1) && (B->col == 1))            )  , "multiply: incorrect length`s", return MAT_FAIL);

    unsigned int r1 = A->row;
    unsigned int r2 = B->row;
    unsigned int c1 = A->col;
    unsigned int c2 = B->col;

    if ((r1 == 1) && (c1 == 1)) {
        scalarmultiply(B, Dest, A->data[0][0]);
        return MAT_SUCC;
    } else if ((r2 == 1) && (c2 == 1)) {
        scalarmultiply(A, Dest, B->data[0][0]);
        return MAT_SUCC;
    }

    register MAT_TYPE sum = (MAT_TYPE)0;
    for (unsigned int i = 0; i < r1; ++i) {
        for (unsigned int j = 0; j < c2; ++j) {
            for (unsigned int k = 0; k < r2; ++k) {
                sum += A->data[i][k] * B->data[k][j];
            }
            Dest->data[i][j] = sum;
            sum = 0;
        }
    }
    return MAT_SUCC;
}

int transpose(Mat* A, Mat* Dest)
{
    M_Assert_Break((!A || !Dest), "transpose: incorrect input values", return MAT_FAIL);
    M_Assert_Break(((A->row > Dest->col)   ||   (A->col > Dest->row)), "transpose: incorrect destination matrix LENGTH", return MAT_FAIL);
    for(unsigned int i = 0; i < A->row; ++i) {
        for(unsigned int j = 0; j < A->col; ++j){
            Dest->data[j][i] = A->data[i][j];
        }
    }
    return MAT_SUCC;
}

Mat* copyValue(Mat* A)
{
    M_Assert_Break((A == NULL), "copyValue: incorrect input", return NULL);
    Mat* B = matrixCreate(A->row,A->col);
    for(unsigned int i = 0; i < A->row; ++i) {
        for(unsigned int j = 0; j < A->col; ++j){
            B->data[i][j] = A->data[i][j];
        }
    }
    return B;
}


int gluInvertMatrix4x4(Mat* A, Mat* Dest) // hardness function: 200
{
    M_Assert_Break((A == NULL || Dest == NULL), "gluInvertMatrix4x4: incorrect input", return MAT_FAIL);
    M_Assert_Break((A->row != 4 || A->col != 4 || Dest->col < 4 || Dest->row < 4), "gluInvertMatrix4x4: incorrect length", return MAT_FAIL);

    float inv[4][4], det;
    inv[0][0]  =  A->data[1][1] * A->data[2][2] * A->data[3][3] - A->data[1][1] * A->data[2][3] * A->data[3][2] - A->data[2][1] * A->data[1][2] * A->data[3][3] + A->data[2][1] * A->data[1][3] * A->data[3][2] + A->data[3][1] * A->data[1][2] * A->data[2][3] - A->data[3][1] * A->data[1][3] * A->data[2][2];
    inv[1][0]  = -A->data[1][0] * A->data[2][2] * A->data[3][3] + A->data[1][0] * A->data[2][3] * A->data[3][2] + A->data[2][0] * A->data[1][2] * A->data[3][3] - A->data[2][0] * A->data[1][3] * A->data[3][2] - A->data[3][0] * A->data[1][2] * A->data[2][3] + A->data[3][0] * A->data[1][3] * A->data[2][2];
    inv[2][0]  =  A->data[1][0] * A->data[2][1]  * A->data[3][3] - A->data[1][0] * A->data[2][3] * A->data[3][1] - A->data[2][0] * A->data[1][1] * A->data[3][3] + A->data[2][0] * A->data[1][3] * A->data[3][1] + A->data[3][0] * A->data[1][1] * A->data[2][3] - A->data[3][0] * A->data[1][3] * A->data[2][1];
    inv[3][0] = -A->data[1][0] * A->data[2][1]  * A->data[3][2] + A->data[1][0] * A->data[2][2] * A->data[3][1] + A->data[2][0] * A->data[1][1] * A->data[3][2] - A->data[2][0] * A->data[1][2] * A->data[3][1] - A->data[3][0] * A->data[1][1] * A->data[2][2] + A->data[3][0] * A->data[1][2] * A->data[2][1];
    inv[0][1]  = -A->data[0][1] * A->data[2][2] * A->data[3][3] + A->data[0][1] * A->data[2][3] * A->data[3][2] + A->data[2][1] * A->data[0][2] * A->data[3][3] - A->data[2][1] * A->data[0][3] * A->data[3][2] - A->data[3][1] * A->data[0][2] * A->data[2][3] + A->data[3][1] * A->data[0][3] * A->data[2][2];
    inv[1][1]  =  A->data[0][0] * A->data[2][2] * A->data[3][3] - A->data[0][0] * A->data[2][3] * A->data[3][2] - A->data[2][0] * A->data[0][2] * A->data[3][3] + A->data[2][0] * A->data[0][3] * A->data[3][2] + A->data[3][0] * A->data[0][2] * A->data[2][3] - A->data[3][0] * A->data[0][3] * A->data[2][2];
    inv[2][1]  = -A->data[0][0] * A->data[2][1]  * A->data[3][3] + A->data[0][0] * A->data[2][3] * A->data[3][1] + A->data[2][0] * A->data[0][1] * A->data[3][3] - A->data[2][0] * A->data[0][3] * A->data[3][1] - A->data[3][0] * A->data[0][1] * A->data[2][3] + A->data[3][0] * A->data[0][3] * A->data[2][1];
    inv[3][1] =  A->data[0][0] * A->data[2][1]  * A->data[3][2] - A->data[0][0] * A->data[2][2] * A->data[3][1] - A->data[2][0] * A->data[0][1] * A->data[3][2] + A->data[2][0] * A->data[0][2] * A->data[3][1] + A->data[3][0] * A->data[0][1] * A->data[2][2] - A->data[3][0] * A->data[0][2] * A->data[2][1];
    inv[0][2]  =  A->data[0][1] * A->data[1][2]  * A->data[3][3] - A->data[0][1] * A->data[1][3]  * A->data[3][2] - A->data[1][1] * A->data[0][2] * A->data[3][3] + A->data[1][1] * A->data[0][3] * A->data[3][2] + A->data[3][1] * A->data[0][2] * A->data[1][3]  - A->data[3][1] * A->data[0][3] * A->data[1][2];
    inv[1][2]  = -A->data[0][0] * A->data[1][2]  * A->data[3][3] + A->data[0][0] * A->data[1][3]  * A->data[3][2] + A->data[1][0] * A->data[0][2] * A->data[3][3] - A->data[1][0] * A->data[0][3] * A->data[3][2] - A->data[3][0] * A->data[0][2] * A->data[1][3]  + A->data[3][0] * A->data[0][3] * A->data[1][2];
    inv[2][2] =  A->data[0][0] * A->data[1][1]  * A->data[3][3] - A->data[0][0] * A->data[1][3]  * A->data[3][1] - A->data[1][0] * A->data[0][1] * A->data[3][3] + A->data[1][0] * A->data[0][3] * A->data[3][1] + A->data[3][0] * A->data[0][1] * A->data[1][3]  - A->data[3][0] * A->data[0][3] * A->data[1][1];
    inv[3][2] = -A->data[0][0] * A->data[1][1]  * A->data[3][2] + A->data[0][0] * A->data[1][2]  * A->data[3][1] + A->data[1][0] * A->data[0][1] * A->data[3][2] - A->data[1][0] * A->data[0][2] * A->data[3][1] - A->data[3][0] * A->data[0][1] * A->data[1][2]  + A->data[3][0] * A->data[0][2] * A->data[1][1];
    inv[0][3]  = -A->data[0][1] * A->data[1][2]  * A->data[2][3] + A->data[0][1] * A->data[1][3]  * A->data[2][2] + A->data[1][1] * A->data[0][2] * A->data[2][3] - A->data[1][1] * A->data[0][3] * A->data[2][2] - A->data[2][1]  * A->data[0][2] * A->data[1][3]  + A->data[2][1]  * A->data[0][3] * A->data[1][2];
    inv[1][3]  =  A->data[0][0] * A->data[1][2]  * A->data[2][3] - A->data[0][0] * A->data[1][3]  * A->data[2][2] - A->data[1][0] * A->data[0][2] * A->data[2][3] + A->data[1][0] * A->data[0][3] * A->data[2][2] + A->data[2][0]  * A->data[0][2] * A->data[1][3]  - A->data[2][0]  * A->data[0][3] * A->data[1][2];
    inv[2][3] = -A->data[0][0] * A->data[1][1]  * A->data[2][3] + A->data[0][0] * A->data[1][3]  * A->data[2][1]  + A->data[1][0] * A->data[0][1] * A->data[2][3] - A->data[1][0] * A->data[0][3] * A->data[2][1]  - A->data[2][0]  * A->data[0][1] * A->data[1][3]  + A->data[2][0]  * A->data[0][3] * A->data[1][1];
    inv[3][3] =  A->data[0][0] * A->data[1][1]  * A->data[2][2] - A->data[0][0] * A->data[1][2]  * A->data[2][1]  - A->data[1][0] * A->data[0][1] * A->data[2][2] + A->data[1][0] * A->data[0][2] * A->data[2][1]  + A->data[2][0]  * A->data[0][1] * A->data[1][2]  - A->data[2][0]  * A->data[0][2] * A->data[1][1];

    det = A->data[0][0] * inv[0][0] + A->data[0][1] * inv[1][0] + A->data[0][2] * inv[2][0] + A->data[0][3] * inv[3][0];
    if (det == 0) {
        return MAT_FAIL;
    }
    det = 1.0f / det;

    for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            Dest->data[i][j] = inv[i][j] * det;
        }
    }

    return MAT_SUCC;
}

int gluInvertMatrix4x4_fastest(Mat* restrict A, Mat* restrict Dest) // destination must not equal to A
{
    M_Assert_Break((A == NULL || Dest == NULL), "gluInvertMatrix4x4_fastest: incorrect input", return MAT_FAIL);
    M_Assert_Break((A->row != 4 || A->col != 4 || Dest->col < 4 || Dest->row < 4), "gluInvertMatrix4x4_fastest: incorrect length", return MAT_FAIL);
    M_Assert_Break((A == Dest), "gluInvertMatrix4x4_fastest: input matrix have equals adresses - undefined behavior", return MAT_FAIL);

    float s[6];
    float c[6];

    s[0] = A->data[0][0]*A->data[1][1] - A->data[1][0]*A->data[0][1];
    s[1] = A->data[0][0]*A->data[1][2] - A->data[1][0]*A->data[0][2];
    s[2] = A->data[0][0]*A->data[1][3] - A->data[1][0]*A->data[0][3];
    s[3] = A->data[0][1]*A->data[1][2] - A->data[1][1]*A->data[0][2];
    s[4] = A->data[0][1]*A->data[1][3] - A->data[1][1]*A->data[0][3];
    s[5] = A->data[0][2]*A->data[1][3] - A->data[1][2]*A->data[0][3];

    c[0] = A->data[2][0]*A->data[3][1] - A->data[3][0]*A->data[2][1];
    c[1] = A->data[2][0]*A->data[3][2] - A->data[3][0]*A->data[2][2];
    c[2] = A->data[2][0]*A->data[3][3] - A->data[3][0]*A->data[2][3];
    c[3] = A->data[2][1]*A->data[3][2] - A->data[3][1]*A->data[2][2];
    c[4] = A->data[2][1]*A->data[3][3] - A->data[3][1]*A->data[2][3];
    c[5] = A->data[2][2]*A->data[3][3] - A->data[3][2]*A->data[2][3];

    /* Assumes it is invertible */
    float idet = 1.0f/( s[0]*c[5]-s[1]*c[4]+s[2]*c[3]+s[3]*c[2]-s[4]*c[1]+s[5]*c[0] );
    if(idet == 0) {
        return MAT_FAIL;
    }

    Dest->data[0][0] = ( A->data[1][1] * c[5] - A->data[1][2] * c[4] + A->data[1][3] * c[3]) * idet;
    Dest->data[0][1] = (-A->data[0][1] * c[5] + A->data[0][2] * c[4] - A->data[0][3] * c[3]) * idet;
    Dest->data[0][2] = ( A->data[3][1] * s[5] - A->data[3][2] * s[4] + A->data[3][3] * s[3]) * idet;
    Dest->data[0][3] = (-A->data[2][1] * s[5] + A->data[2][2] * s[4] - A->data[2][3] * s[3]) * idet;

    Dest->data[1][0] = (-A->data[1][0] * c[5] + A->data[1][2] * c[2] - A->data[1][3] * c[1]) * idet;
    Dest->data[1][1] = ( A->data[0][0] * c[5] - A->data[0][2] * c[2] + A->data[0][3] * c[1]) * idet;
    Dest->data[1][2] = (-A->data[3][0] * s[5] + A->data[3][2] * s[2] - A->data[3][3] * s[1]) * idet;
    Dest->data[1][3] = ( A->data[2][0] * s[5] - A->data[2][2] * s[2] + A->data[2][3] * s[1]) * idet;

    Dest->data[2][0] = ( A->data[1][0] * c[4] - A->data[1][1] * c[2] + A->data[1][3] * c[0]) * idet;
    Dest->data[2][1] = (-A->data[0][0] * c[4] + A->data[0][1] * c[2] - A->data[0][3] * c[0]) * idet;
    Dest->data[2][2] = ( A->data[3][0] * s[4] - A->data[3][1] * s[2] + A->data[3][3] * s[0]) * idet;
    Dest->data[2][3] = (-A->data[2][0] * s[4] + A->data[2][1] * s[2] - A->data[2][3] * s[0]) * idet;

    Dest->data[3][0] = (-A->data[1][0] * c[3] + A->data[1][1] * c[1] - A->data[1][2] * c[0]) * idet;
    Dest->data[3][1] = ( A->data[0][0] * c[3] - A->data[0][1] * c[1] + A->data[0][2] * c[0]) * idet;
    Dest->data[3][2] = (-A->data[3][0] * s[3] + A->data[3][1] * s[1] - A->data[3][2] * s[0]) * idet;
    Dest->data[3][3] = ( A->data[2][0] * s[3] - A->data[2][1] * s[1] + A->data[2][2] * s[0]) * idet;

    return MAT_SUCC;
}
