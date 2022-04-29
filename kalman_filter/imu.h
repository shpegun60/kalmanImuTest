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
    float grav[3];              // result calculated gravity vector
    float const_u;              // angle error filtration constant coefficient u (set once at the beginning) Must be: in conditions of static motion >= 0.9, in conditions of dynamic motion <= 0.001
    float calibrationGrav[3];   // calibration gravity constant on all axes (set once at the beginning)

    // temp data
    float angleErr;             // error angle between calculated accelerometer data and measured accelerometer
    float rotateAxisErr[3];     // error axis between calculated accelerometer data and measured accelerometer
    float grav_norm[3];         // normalized gravity vector
    Quaternion q_ae;            // error rotation quaternion
    Quaternion q_a;             // result estimated quaternion from kalman

    // other temp data for help`s calculation
    float Tx;
    float Ty;
    float Tz;
    float halfT;
    float recipNorm;
} IMU10Dof;

typedef struct
{
    float dt;

    float g[3]; // [ radians/sec ]   g[0] ==> x; g[1] ==> y; g[2] ==> z;
    float a[3]; // [ m / sec^2 ]  a[0] -> x; a[1] -> y; a[2] -> z; earth gravity ==> z
    float m[3]; // [ones Gauss]    m[0] -> x; m[1] -> y; m[2] -> z; (on this time not used)

} IMUinput;


IMU10Dof* imuCreate(float const_u, float* gravityConst, MAT_TYPE dt, MAT_TYPE* Q10x10, MAT_TYPE* R4x4);
int imuProceed(IMU10Dof* imu, IMUinput * data);



#endif /* __IMU */

