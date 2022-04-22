#include <iostream>

extern "C" {
#include "kalman_filter.h"
#include "test_quaternion.h"

#include "matrix.h"
#include "matrix_test.h"
}

#include <time.h>

using namespace std;

void kalmanInit(Mat* X_est, Mat* P_est, Mat* F, Mat* F_t, Mat* G, Mat* Q, Mat* R, Mat* H, Mat* H_t)
{
    return;
}

int main()
{
    test_quaternion();

    testMatrix();

    KalmanFilter* kalman = kalmanCreate(kalmanInit, 10, 4, 3);


    for(int i =0; i < 100; ++i) {
        kalmanPredict(kalman);
        kalmanUpdate(kalman);
    }

        printf("---------------- KALMAN PRINT----------------------------------- \n\n");
        printkalman(kalman);
        printf("\nsizeof Kalman: %d\n", (int)sizeof(KalmanFilter));

    return 0;
}
