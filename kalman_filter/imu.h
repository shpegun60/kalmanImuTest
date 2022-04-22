/**
 * @file    IMU.h
 * @brief
 * @date
 */

#ifndef __IMU
#define __IMU

#include "kalman_filter.h"

typedef struct
{
    KalmanFilter* kalman;
    float grav[3];
    float Tx;
    float Ty;
    float Tz;
    float halfT;
    float quadHalfT;
} IMU10Dof;

IMU10Dof* imuCreate(MAT_TYPE dt, MAT_TYPE* Q10x10, MAT_TYPE* R4x4);

typedef struct
{
    float dt;

    float gx; // rad/s
    float gy; // rad/s
    float gz; // rad/s

    float ax;
    float ay;
    float az; // normal gravity component

    float mx;
    float my;
    float mz;

} IMUinput;



#endif /* __IMU */

