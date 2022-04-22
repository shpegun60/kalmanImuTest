#include "matrix_test.h"
#include "matrix.h"
#include <stdio.h>

void testMatrix(void) {
    /* first step
    The original of matrix:
    [	1.00000		0.00000		4.00000		-6.00000
        2.00000		5.00000		0.00000		3.00000
        -1.00000	2.00000		3.00000		5.00000
        2.00000		1.00000		-2.00000	3.00000	]
     */

    Mat * INV = matrixCreate(4, 4);
    Mat * INV_dest = matrixCreate(4, 4);
    INV->data[0][0] = 1;
    INV->data[0][1] = 0;
    INV->data[0][2] = 4;
    INV->data[0][3] = -6;

    INV->data[1][0] = 2;
    INV->data[1][1] = 5;
    INV->data[1][2] = 0;
    INV->data[1][3] = 3;

    INV->data[2][0] = -1;
    INV->data[2][1] = 2;
    INV->data[2][2] = 3;
    INV->data[2][3] = 5;

    INV->data[3][0] = 2;
    INV->data[3][1] = 1;
    INV->data[3][2] = -2;
    INV->data[3][3] = 3;

    showmat(INV, (char *)"The original of matrix:");
    int boolr = gluInvertMatrix4x4_fastest(INV, INV_dest);
    showmat(INV_dest, (char *)"\nThe inverse of matrix:");
    printf("\n bool = %d \n\n", boolr);

    /* true result
    The inverse of matrix:
    [	0.23270		-0.11950	0.03774		0.52201
        -0.08176	0.29874		-0.09434	-0.30503
        0.16352		-0.09748	0.18868		0.11006
        -0.01887	-0.08491	0.13208		0.16038	]
    */

    // second step-------------------------------------
    Mat * A = matrixCreate(3, 2);
    Mat * B = matrixCreate(2, 3);
    Mat * D = matrixCreate(1, 1);
    Mat * C = createResultMulMatrix(A, B);
    Mat * T = createResultTransMatrix(B);

    //Mat * Drr = matrixCreate(0, 0);

    A->data[0][0] = 2;
    A->data[0][1] = 3;
    A->data[1][0] = -5;
    A->data[1][1] = 6;
    A->data[2][0] = 9;
    A->data[2][1] = -7;

    B->data[0][0] = 1;
    B->data[0][1] = -2;
    B->data[0][2] = 0;
    B->data[1][0] = 3;
    B->data[1][1] = 4;
    B->data[1][2] = -5;

    D->data[0][0] = 2;

    transpose(B, T);


    multiply(A, B, C);

    showmat(A, (char *)"A:");
    showmat(B, (char *)"B:");
    showmat(C, (char *)"A*B = C:");
    showmat(T, (char *)"B T:");


    Mat* Copy = copyValue(A);
    showmat(Copy, (char *)"Copy A:");


    multiply(A, D, A);
    showmat(D, (char *)"D:");
    showmat(A, (char *)"A*D = A:");


    sub(A, T, A);
    showmat(A, (char *)"A-BT = A:");

    destroy_matrix(&A);
    showmat(A, (char *)"destroied matrix A:");

    //transpose(NULL, NULL); // check error message

    //    //-------------------------------------------------
}
