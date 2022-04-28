#include "kalmanImuTest.h"
#include <qdebug.h>

KalmanIMUTest::KalmanIMUTest()
{
    test_quaternion();
    //testMatrix();

    float sR = 0.0015;
    float sQR = 1e-6;
    float sQx = 0.001;
    float sQy = 0.001;
    float sQz = 0.001;

    float dt = 1/1000;

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

    kalman = new KalmanIMU(0.9, 0, Q_init[0], R_init[0]);
}



void KalmanIMUTest::testQuaternionKalman(QCPGraph* graphX, QCPGraph* graphY, QCPGraph* graphZ, QString gyroFileName, QString accFileName)
{

    graphX->setName("alpha");
    graphY->setName("beta");
    graphZ->setName("gamma");

    QFile gyroFile(gyroFileName);
    QFile accFile(accFileName);

    if (!gyroFile.open(QIODevice::ReadOnly)) {
        qDebug() <<"gyro file not open: "<< gyroFile.errorString();
        return;
    }

    if (!accFile.open(QIODevice::ReadOnly)) {
        qDebug() <<endl<<"acc file not open: "<< accFile.errorString();
        return;
    }

    bool ok;
    float euler[3] = {0.0, };
    quint64 start_time = 0;
    quint64 last_time = 0;
    bool timefixed = false;
    bool headPass = false;


    while (!gyroFile.atEnd() && !accFile.atEnd()) {
        QByteArray accLine = accFile.readLine();
        QByteArray gyroLine = gyroFile.readLine();

        if(!headPass) {
            headPass = true;
            continue;
        }

        auto &&accList = accLine.split(',');
        auto &&gyroList = gyroLine.split(',');

        if(!timefixed) {
            start_time = gyroList[0].toULongLong(&ok);
            timefixed = true;
            continue;
        }

        auto t = gyroList[0].toULongLong(&ok) - start_time;
        auto dt = t - last_time;
        last_time = t;

        if(dt == 0) {
            continue;
        }

        auto ax = accList[1].toDouble(&ok);
        auto ay = accList[2].toDouble(&ok);
        auto az = accList[3].toDouble(&ok);

        auto gx = gyroList[1].toDouble(&ok);
        auto gy = gyroList[2].toDouble(&ok);
        auto gz = gyroList[3].toDouble(&ok);

        kalman->KalmanIMUProceed(dt / 1000.0, ax, ay, az, gx, gy, gz, 0, 0, 0);
        Quaternion_toEulerZYX(kalman->getQuaternion(), euler);
        //ui->qplot1->graph(0)->addData(i, y0[i]);

        //kalman->printKalmanTop();
        graphX->addData(t, euler[0] * (180.0f / M_PI));
        graphY->addData(t, euler[1] * (180.0f / M_PI));
        graphZ->addData(t, euler[2] * (180.0f / M_PI));
    }

    graphX->parentPlot()->replot();
    graphY->parentPlot()->replot();
    graphZ->parentPlot()->replot();

}



void KalmanIMUTest::csvXYZToGraph(QString fileName, QCPGraph* graphX, QCPGraph* graphY, QCPGraph* graphZ, quint32 counterStart, double multiplicationConst)
{

    graphX->setName("X");
    graphY->setName("Y");
    graphZ->setName("Z");

    QFile file(fileName);

    if (!file.open(QIODevice::ReadOnly)) {
        qDebug()<<endl<<fileName<<"file not open: "<< file.errorString();
        return;
    }

    bool ok;
    quint64 start_time = 0;
    quint64 last_time = 0;
    bool headPass = false;
    quint32 counter = counterStart;
    double startX = 0.0f;
    double startY = 0.0f;
    double startZ = 0.0f;


    while (!file.atEnd()) {
        QByteArray line = file.readLine();

        if(!headPass) {
            headPass = true;
            continue;
        }

        auto &&list = line.split(',');

        if(counter) {
            start_time = list[0].toULongLong(&ok);

            if(counterStart > 1) {
                startX += list[1].toDouble(&ok);
                startY += list[2].toDouble(&ok);
                startZ += list[3].toDouble(&ok);
                --counter;

                if(counter == 0) {
                    startX /=(float)counterStart;
                    startY /=(float)counterStart;
                    startZ /=(float)counterStart;
                }
            } else {
                startX = 0.0f;
                startY = 0.0f;
                startZ = 0.0f;
                --counter;
            }
            continue;
        }

        auto t = list[0].toULongLong(&ok) - start_time;
        auto dt = t - last_time;
        last_time = t;

        if(dt == 0) {
            continue;
        }

        auto dx = (list[1].toDouble(&ok) - startX)  * multiplicationConst;
        auto dy = (list[2].toDouble(&ok) - startY)  * multiplicationConst;
        auto dz = (list[3].toDouble(&ok) - startZ)  * multiplicationConst;


        graphX->addData(t, dx);
        graphY->addData(t, dy);
        graphZ->addData(t, dz);
    }

    graphX ->parentPlot()->replot();
    graphY ->parentPlot()->replot();
    graphZ ->parentPlot()->replot();
}
