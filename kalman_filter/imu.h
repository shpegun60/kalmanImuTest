/**
 * @file    IMU.h
 * @brief
 * @date
 */

#ifndef __IMU
#define __IMU

#include "kalman_filter.h"
#include "quaternion.h"

typedef struct
{
    KalmanFilter* kalman;
    float grav[3];
    float const_u;

    float angleErr;
    float rotateAxisErr[3];
    float grav_norm[3];
    Quaternion q_ae;
    Quaternion q_a;

    float Tx;
    float Ty;
    float Tz;
    float halfT;
    float quadHalfT;
    float recipNorm;
} IMU10Dof;

typedef struct
{
    float dt;

    float g[3]; // [radians/sec]   g[0] ==> x; g[1] ==> y; g[2] ==> z;
    float a[3]; // [ones Gravity]  a[0] -> x; a[1] -> y; a[2] -> z; earth gravity ==> z
    float m[3]; // [ones Gauss]    m[0] -> x; m[1] -> y; m[2] -> z;

} IMUinput;


IMU10Dof* imuCreate(float const_u, MAT_TYPE dt, MAT_TYPE* Q10x10, MAT_TYPE* R4x4);
int imuProceed(IMU10Dof* imu, IMUinput * data);



#endif /* __IMU */

