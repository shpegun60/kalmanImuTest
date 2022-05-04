#ifndef IMU_TEST_H
#define IMU_TEST_H

#include "kalmanImuTop.h"
#include <qcustomplot.h>
#include <qstring.h>


typedef enum {
    ANGLE_DATA,
    GRAVITY_DATA,
    LINEAR_ACCERARATION_DATA
} RequestTests;



class KalmanIMUTest
{
public:
    KalmanIMUTest();

    void testQuaternionKalman(QCPGraph* graphX, QCPGraph* graphY, QCPGraph* graphZ, QString gyroFileName, QString accFileName, RequestTests testType);
    void csvXYZToGraph(QString fileName, QCPGraph* graphX,QCPGraph* graphY, QCPGraph* graphZ, quint32 counterStart = 0, double multiplicationConst = 1.0f);

    void calibrateAccelGyroFromFile(QCPGraph *accX, QCPGraph *accY, QCPGraph *accZ, QCPGraph *gyroX, QCPGraph *gyroY, QCPGraph *gyroZ, QString gyroFileName, QString accFileName, double gravityConstant, int meanIterations, double acel_deadzone, double gyro_deadzone);

private:
    void calculateMeans(double ax, double ay, double az, double gx, double gy, double gz, int meanIterations);

private:
    KalmanIMU* kalman;


    // calibration values
    double buff_ax = 0, buff_ay = 0, buff_az = 0, buff_gx = 0, buff_gy = 0, buff_gz = 0;
    double mean_ax = 0, mean_ay = 0, mean_az = 0, mean_gx = 0, mean_gy = 0, mean_gz = 0;
    double ax_offset = 0, ay_offset = 0, az_offset = 0, gx_offset = 0, gy_offset = 0, gz_offset = 0;
    int iterator = 0;
    int state = 0;
    int NOC = 0;

    double mean_ax_best = 0, mean_ay_best = 0, mean_az_best = 0, mean_gx_best = 0, mean_gy_best = 0, mean_gz_best = 0;
    double ax_offset_best = 0, ay_offset_best = 0, az_offset_best = 0, gx_offset_best = 0, gy_offset_best = 0, gz_offset_best = 0;
    double bestSum = 0.0;

};

#endif // IMU_TEST_H
