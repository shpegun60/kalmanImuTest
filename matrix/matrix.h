#ifndef __MATRIX
#define __MATRIX


#define MAT_SUCC 1
#define MAT_FAIL -1

#ifndef NULL
    #define NULL ((void *)0)
#endif

/* current representation of a matrix in my mind  */

#define MAT_TYPE float

typedef struct {
    unsigned int row;
    unsigned int col;
    MAT_TYPE **data;// data[row][col]
} Mat;

int showmat(Mat* A, char * name);

/* make a zero matrix of given dimensions */
Mat *matrixCreate(unsigned int r, unsigned int c);
int matrixInitFromArr(Mat* A, MAT_TYPE* arr);

int matrixInitFromArr_T(Mat* A, MAT_TYPE* arr);
int getResultRawSize(Mat* A, Mat* B, unsigned int * resRaw, unsigned int * resCall);
Mat *createResultMulMatrix(Mat* A, Mat* B);
Mat *createResultTransMatrix(Mat* A);

/* free memory associated with the matrix  */
int destroy_matrix(Mat **m);

/* enter 1s along the main diagonal */
Mat *eye(unsigned int length);
Mat* ones(unsigned int r, unsigned int c, MAT_TYPE d);
int add(Mat* A, Mat* B, Mat* Dest);

int sub(Mat* A, Mat* B, Mat* Dest);
int scalarmultiply(Mat* A, Mat* Dest, MAT_TYPE scalar);

int multiply(Mat* A, Mat* B, Mat* Dest);

int transpose(Mat* A, Mat* Dest);
Mat* copyValue(Mat* A);
int matrixCopy(Mat* A, Mat* Dest);

int gluInvertMatrix4x4(Mat* A, Mat* Dest);
int gluInvertMatrix4x4_fastest(Mat* A, Mat* Dest); // destination must not equal to A


#endif // __MATRIX
