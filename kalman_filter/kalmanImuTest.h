#ifndef IMU_TEST_H
#define IMU_TEST_H

#include "kalmanImuTop.h"
#include <qcustomplot.h>
#include <qstring.h>

class KalmanIMUTest
{
public:
    KalmanIMUTest();

    void testQuaternionKalman(QCPGraph* graphX, QCPGraph* graphY, QCPGraph* graphZ, QString gyroFileName, QString accFileName);
    void csvXYZToGraph(QString fileName, QCPGraph* graphX,QCPGraph* graphY, QCPGraph* graphZ, quint32 counterStart = 0, double multiplicationConst = 1.0f);

private:
    KalmanIMU* kalman;
};

#endif // IMU_TEST_H
