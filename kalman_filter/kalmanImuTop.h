#ifndef IMU_H
#define IMU_H

extern "C" {
#include "imu.h"
#include "test_quaternion.h"
#include "matrix_test.h"
}

#include <qstring.h>

class KalmanIMU
{
public:
    KalmanIMU(float const_u, float dt, float* Q10x10, float* R4x4);
    int KalmanIMUProceed(IMUinput * data);
    int KalmanIMUProceed(double dt, double ax, double ay, double az, double gx, double gy, double gz, double mx, double my, double mz);


    Mat * getResultData();
    float *getGravityData();
    Quaternion *getQuaternion();

    void showMatrixTop(Mat *, QString name);
    void printKalmanTop();

private:
    IMU10Dof* imu;
};

#endif // IMU_H