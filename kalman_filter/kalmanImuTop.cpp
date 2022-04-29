#include "kalmanImuTop.h"
#include <qdebug.h>

KalmanIMU::KalmanIMU(float const_u, float* gravityConst, float dt, float *Q10x10, float *R4x4)
{
    imu = imuCreate(const_u, gravityConst, dt, Q10x10, R4x4);
    //    /printkalman(imu->kalman);
}

int KalmanIMU::KalmanIMUProceed(IMUinput *data)
{
    return imuProceed(imu, data);
}

int KalmanIMU::KalmanIMUProceed(double dt, double ax, double ay, double az, double gx, double gy, double gz, double mx, double my, double mz)
{
    IMUinput data;

    data.dt = dt;
    data.a[0] = ax;
    data.a[1] = ay;
    data.a[2] = az;

    data.g[0] = gx;
    data.g[1] = gy;
    data.g[2] = gz;

    data.m[0] = mx;
    data.m[1] = my;
    data.m[2] = mz;

//    char txt[500];
//    sprintf(txt, "\nNEW ITERATION gx: %0.5f  gy: %0.5f gz: %0.5f  ax: %0.5f  ay: %0.5f  az: %0.5f  T: %0.5f -------------------------------------------------", gx,gy,gz,ax,ay,az,dt);
//    printf("%s", (char*)txt);
//    sprintf(txt, "\nGravity[X]: %0.5f Gravity[Y]: %0.5f Gravity[Z]: %0.5f \n", imu->grav[0], imu->grav[1],imu->grav[1]);
//    printf("%s", (char*)txt);

    return imuProceed(imu, &data);
}

Mat *KalmanIMU::getResultData()
{
    return imu->kalman->X_est;
}

float *KalmanIMU::getGravityData()
{
    return imu->grav;
}

Quaternion *KalmanIMU::getQuaternion()
{
    return &imu->q_a;
}

float *KalmanIMU::getGravity()
{
    return imu->grav;
}


void KalmanIMU::printKalmanTop()
{
    printkalman(imu->kalman);
}


