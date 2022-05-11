/**
 * @file    IMU.h
 * @brief
 * @date
 */

#ifndef __IMU
#define __IMU

#include "kalman_filter.h"
#include "quaternion.h"

/*
There are two main global reference frames based on the local tangent plane:
    -NED defines the X-, Y-, and Z-axis colinear to the geographical North, East, and Down directions, respectively.
    -ENU defines the X-, Y-, and Z-axis colinear to the geographical East, North, and Up directions, respectively.

    if NED gravity vector G = [0, 0, -1]
    if ENU gravity vector G = [0, 0,  1]


    if NED geomagnetic field vector R = [??????]
    if ENU geomagnetic field vector R = [??????]
*/

typedef enum {
    NED,
    ENU
} GlobalCoordinate;

typedef struct {
    float dt_init_sec;
    float accConst_u;               // accelerometer influence factor to gyro quaternion
    float magConst_u;               // magnetometer influence factor to gyro quaternion
    float gravityConstVect[3];      // gravity field for each axes [x, y, z] (must be positive values)==> [9.81, 9.81, 9.81] = ABS(GRAVITY)

    // biases
    float accelBiasVect[3];         // accelerometer bias vector [x, y, z] ==> (accel + bias)
    float gyroBiasVect[3];          // gyroscope bias vector [x, y, z] ==> (gyro + bias)

    // variances
    float gyroVarianceVect[3];   // gyroscope variance vector [x, y, z]
    float accVarianceVect[3];    // accelerometer variance vector [x, y, z]
    float magVarianceVect[3];    // magnetometer variance vector [x, y, z]

    // coordinate system
    GlobalCoordinate coordinateType; // coordinate system
} IMUInit_struct;

typedef struct
{
    float dt_sec;
    float g[3]; // [ radians/sec ]   g[0] ==> x; g[1] ==> y; g[2] ==> z;
    float a[3]; // [ m / sec^2 ]  a[0] -> x; a[1] -> y; a[2] -> z; earth gravity ==> z
    float m[3]; // [ones Gauss]    m[0] -> x; m[1] -> y; m[2] -> z; (on this time not used)
} IMUinput;


typedef struct IMU IMU10Dof;
struct IMU
{
    KalmanFilter* kalman;
    float grav[3];              // result calculated gravity vector
    float accConst_u;           // angle error filtration constant coefficient u for accelerometer (set once at the beginning) Must be: in conditions of static motion >= 0.9, in conditions of dynamic motion <= 0.001
    float calibrationGrav[3];   // calibration gravity constant on all axes (set once at the beginning)
    float accelBiasVect[3];   // acceleration bias vector
    float gyroBiasVect[3];   // gyroscope bias vector

    // temp data
    float rotateAxisErr[3];     // error axis between calculated accelerometer data and measured accelerometer
    float grav_norm[3];         // normalized gravity vector
    Quaternion q_ae;            // error rotation quaternion
    Quaternion q_a;             // result estimated quaternion from kalman
    Quaternion q_i;             // identity quaternion for filtration

    // other temp data for help`s calculation
    float Tx;
    float Ty;
    float Tz;
    float halfT;
    float recipNorm;

    // to update Q_k------------------------------------------------------
    Mat* Q_quat;
    Mat* G;
    Mat* G_t;
    Mat* Noise;
    Mat* NOISE_RES;
    float constT_4;

    // to update R --------------------------------------------------------
    Mat* J;
    Mat* J_t;
    Mat* Noise_acc;
    Mat* NOISE_R_RES;


    void (*accDeltaQuaterFinder) (IMU10Dof* imu, IMUinput * data); // accelerometer delta quaternion finder for choosed coordinates

};



IMU10Dof* imuCreate(IMUInit_struct* init);
int imuProceed(IMU10Dof* imu, IMUinput* data);


void accDeltaQuaterFinder_NED(IMU10Dof* imu, IMUinput * data);
void accDeltaQuaterFinder_ENU(IMU10Dof* imu, IMUinput * data);

#endif /* __IMU */

