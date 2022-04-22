#include <iostream>

extern "C" {
#include "kalman_filter.h"
#include "test_quaternion.h"

#include "matrix.h"
#include "matrix_test.h"
#include "imu.h"
}

#include <time.h>

using namespace std;

int main()
{
    test_quaternion();

    testMatrix();

    MAT_TYPE Q_init[10][10] =
    {
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  0}
    };

    MAT_TYPE R_init[4][4] =
    {
        {0, 0,  0,  0},
        {0, 0,  0,  0},
        {0, 0,  0,  0},
        {0, 0,  0,  0}
    };

    IMU10Dof* imu = imuCreate(10, Q_init[0], R_init[0]);



    printf("---------------- KALMAN PRINT----------------------------------- \n\n");
    printkalman(imu->kalman);
    printf("\nsizeof Kalman: %d\n", (int)sizeof(KalmanFilter));

    return 0;
}
