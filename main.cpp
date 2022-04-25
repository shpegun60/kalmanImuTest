#include <iostream>
#include <math.h>

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

    float sR = 0.0015;
    float sQR = 1;
    float sQx = 0.001;
    float sQy = 0.001;
    float sQz = 0.001;

    float dt = 10;

    float T11 = (dt * dt * dt * dt) / 4.0;
    float T12 = (dt * dt * dt) / 2.0;
    float T21 = (dt * dt * dt) / 2.0;
    float T22 = (dt * dt);

    MAT_TYPE Q_init[10][10] =
    {
        {T11 * sQx,     T12 * sQx,              0,                      0,                              0,                          0,                  0,  0,  0,  0},
        {T21 * sQx,     T22 * sQx,              0,                      0,                              0,                          0,                  0,  0,  0,  0},
        {0,             0,                      T11 * sQy,              T12 * sQx,                      0,                          0,                  0,  0,  0,  0},
        {0,             0,                      T21 * sQy,              T22 * sQy,                      0,                          0,                  0,  0,  0,  0},
        {0,             0,                      0,                      0,                              T11 * sQz,                  T12 * sQz,          0,  0,  0,  0},
        {0,             0,                      0,                      0,                              T21 * sQz,                  T22 * sQz,          0,  0,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  sQR,  0,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  sQR,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  0,  sQR,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  0,  0,  sQR}
    };

    MAT_TYPE R_init[4][4] =
    {
        {sR, 0,  0,  0},
        {0, sR,  0,  0},
        {0, 0,  sR,  0},
        {0, 0,  0,  sR}
    };

    IMU10Dof* imu = imuCreate(0.9, dt, Q_init[0], R_init[0]);



    printf("---------------- KALMAN PRINT----------------------------------- \n\n");
    printkalman(imu->kalman);
    printf("\nsizeof Kalman: %d\n", (int)sizeof(KalmanFilter));

    return 0;
}
