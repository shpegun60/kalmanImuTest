#include "kalmanImuTest.h"
#include <qdebug.h>


#define FILTER_ALPHA 0.95f

KalmanIMUTest::KalmanIMUTest()
{
    test_quaternion();
    //testMatrix();

    IMUInit_struct kalmanInit;
    kalmanInit.dt_init_sec = 0.0;
    kalmanInit.accConst_u = 1.0;
    kalmanInit.magConst_u = 1.0;

    kalmanInit.gravityConstVect[0] = 9.81;
    kalmanInit.gravityConstVect[1] = 9.81;
    kalmanInit.gravityConstVect[2] = 9.81;

//    kalmanInit.accelBiasVect[0] = 2.63679;
//    kalmanInit.accelBiasVect[1] = -0.064461;
//    kalmanInit.accelBiasVect[2] =  -0.330342;

//    kalmanInit.gyroBiasVect[0] = -0.000133188;
//    kalmanInit.gyroBiasVect[1] = 2.54972e-05;
//    kalmanInit.gyroBiasVect[2] = -9.72103e-06;

    kalmanInit.accelBiasVect[0] = 0.0;
    kalmanInit.accelBiasVect[1] = 0.0;
    kalmanInit.accelBiasVect[2] = 0.0;

    kalmanInit.gyroBiasVect[0] = 0.0;
    kalmanInit.gyroBiasVect[1] = 0.0;
    kalmanInit.gyroBiasVect[2] = 0.0;


    kalmanInit.gyroVarianceVect[0] = 0.0001;
    kalmanInit.gyroVarianceVect[1] = 0.0001;
    kalmanInit.gyroVarianceVect[2] = 0.0001;

    kalmanInit.accVarianceVect[0] = 0.0001;
    kalmanInit.accVarianceVect[1] = 0.0001;
    kalmanInit.accVarianceVect[2] = 0.0001;

    kalmanInit.magVarianceVect[0] = 0.0001;
    kalmanInit.magVarianceVect[1] = 0.0001;
    kalmanInit.magVarianceVect[2] = 0.0001;

    kalmanInit.coordinateType = NED;//ENU;


    kalman = new KalmanIMU(&kalmanInit);
}



void KalmanIMUTest::testQuaternionKalman(QCPGraph* graphX, QCPGraph* graphY, QCPGraph* graphZ, QString gyroFileName, QString accFileName, QString magFileName, RequestTests testType)
{

    graphX->setName("alpha");
    graphY->setName("beta");
    graphZ->setName("gamma");

    QFile gyroFile(gyroFileName);
    QFile accFile(accFileName);
    QFile magFile(magFileName);

    if (!gyroFile.open(QIODevice::ReadOnly)) {
        qDebug() <<"gyro file not open: "<< gyroFile.errorString();
        return;
    }

    if (!accFile.open(QIODevice::ReadOnly)) {
        qDebug() <<endl<<"acc file not open: "<< accFile.errorString();
        return;
    }

    if (!magFile.open(QIODevice::ReadOnly)) {
        qDebug() <<endl<<"mag file not open: "<< magFile.errorString();
        return;
    }

    bool ok;
    float euler[3] = {0.0, };
    quint64 start_time = 0;
    quint64 last_time = 0;
    bool timefixed = false;
    bool headPass = false;

    int b = 0;

    float last_Alpha[3] = {0, 0, 0};
    float d_alpha[3] = {0, 0, 0};

    while (!gyroFile.atEnd() && !accFile.atEnd() && !magFile.atEnd()) {
        QByteArray accLine = accFile.readLine();
        QByteArray gyroLine = gyroFile.readLine();
        QByteArray magLine = magFile.readLine();

        if(!headPass) {
            headPass = true;
            continue;
        }

        auto &&accList = accLine.split(',');
        auto &&gyroList = gyroLine.split(',');
        auto &&magList = magLine.split(',');

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

        auto mx = magList[1].toDouble(&ok);
        auto my = magList[2].toDouble(&ok);
        auto mz = magList[3].toDouble(&ok);

        kalman->KalmanIMUProceed(dt / 1000.0, ax, ay, az, gx, gy, gz, mx, my, mz);


        switch (testType) {

        case(ANGLE_DATA):
            Quaternion_toEulerZYX(kalman->getQuaternion(), euler);

            //kalman->printKalmanTop();
            graphX->addData(t, euler[0] * (180.0f / M_PI));
            graphY->addData(t, euler[1] * (180.0f / M_PI));
            graphZ->addData(t, euler[2] * (180.0f / M_PI));
            break;

        case(ANGLE_VELOCITY_DATA):
            Quaternion_toEulerZYX(kalman->getQuaternion(), euler);

            d_alpha[0] = (1.0f - FILTER_ALPHA) * (((euler[0] - last_Alpha[0]) * (180.0f / M_PI)) / (dt / 1000.0)) + FILTER_ALPHA * d_alpha[0];
            d_alpha[1] = (1.0f - FILTER_ALPHA) * (((euler[1] - last_Alpha[1]) * (180.0f / M_PI)) / (dt / 1000.0)) + FILTER_ALPHA * d_alpha[1];
            d_alpha[2] = (1.0f - FILTER_ALPHA) * (((euler[2] - last_Alpha[2]) * (180.0f / M_PI)) / (dt / 1000.0)) + FILTER_ALPHA * d_alpha[2];


            graphX->addData(t, d_alpha[0]);
            graphY->addData(t, d_alpha[1]);
            graphZ->addData(t, d_alpha[2]);

            last_Alpha[0] = euler[0];
            last_Alpha[1] = euler[1];
            last_Alpha[2] = euler[2];

            break;

        case(GRAVITY_DATA): {
            float * gravity = kalman->getGravity();

            graphX->addData(t, gravity[0]);
            graphY->addData(t, gravity[1]);
            graphZ->addData(t, gravity[2]);
            break;
        }

        case(LINEAR_ACCERARATION_DATA): {
            Mat * linearAcceleration = kalman->getLinearAcceleration();

            graphX->addData(t, linearAcceleration->data[0][0]);
            graphY->addData(t, linearAcceleration->data[1][0]);
            graphZ->addData(t, linearAcceleration->data[2][0]);
            break;


        }

        default:
            break;
        }

        //        b++;
        //        if(b > 1000 ) {
        //            return;
        //        }
    }


    (void)b;

    graphX->parentPlot()->replot();
    graphY->parentPlot()->replot();
    graphZ->parentPlot()->replot();

    gyroFile.close();
    accFile.close();
    magFile.close();

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

    file.close();
}



void KalmanIMUTest::calibrateAccelGyroFromFile(QCPGraph *accX, QCPGraph *accY, QCPGraph *accZ, QCPGraph *gyroX, QCPGraph *gyroY, QCPGraph *gyroZ, QString gyroFileName, QString accFileName, double gravityConstant, int meanIterations, double acel_deadzone, double gyro_deadzone)
{
    accX->setName("accX");
    accY->setName("accY");
    accZ->setName("accZ");

    gyroX->setName("gyroX");
    gyroY->setName("gyroY");
    gyroZ->setName("gyroZ");


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
    bool headPass = false;
    bool timefixed = false;
    quint64 start_time = 0;
    quint64 last_time = 0;



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

        auto ax = accList[1].toDouble(&ok) + ax_offset;
        auto ay = accList[2].toDouble(&ok) + ay_offset;
        auto az = accList[3].toDouble(&ok) + az_offset;

        auto gx = gyroList[1].toDouble(&ok) + gx_offset;
        auto gy = gyroList[2].toDouble(&ok) + gy_offset;
        auto gz = gyroList[3].toDouble(&ok) + gz_offset;

        switch (state) {

        case(0):
            buff_ax = 0.0f; buff_ay = 0.0f; buff_az = 0.0f; buff_gx = 0.0f; buff_gy = 0.0f; buff_gz = 0.0f;
            mean_ax = 0.0f; mean_ay = 0.0f; mean_az = 0.0f; mean_gx = 0.0f; mean_gy = 0.0f; mean_gz = 0.0f;
            ax_offset = 0.0f;  ay_offset = 0.0f; az_offset = 0.0f; gx_offset = 0.0f; gy_offset = 0.0f; gz_offset = 0.0f;
            iterator = 0;
            ++state;
            break;

        case(1):
            calculateMeans(ax, ay, az, gx, gy, gz, meanIterations);
            break;

        case(2):
            ax_offset = -mean_ax / 8.0f;
            ay_offset = -mean_ay / 8.0f;
            az_offset = (gravityConstant - mean_az) / 8.0f;
            gx_offset = -mean_gx / 4.0f;
            gy_offset = -mean_gy / 4.0f;
            gz_offset = -mean_gz / 4.0f;

            mean_ax_best = 0.0; mean_ay_best = 0.0; mean_az_best = 0.0; mean_gx_best = 0.0; mean_gy_best = 0.0; mean_gz_best = 0.0;
            ax_offset_best = 0.0; ay_offset_best = 0.0; az_offset_best = 0.0; gx_offset_best = 0.0; gy_offset_best = 0.0; gz_offset_best = 0.0;
            bestSum = abs(mean_ax) + abs(mean_ay) + abs(gravityConstant - mean_az) + abs(mean_gx) + abs(mean_gy) + abs(mean_gz);

            buff_ax = 0.0f; buff_ay = 0.0f; buff_az = 0.0f; buff_gx = 0.0f; buff_gy = 0.0f; buff_gz = 0.0f;
            iterator = 0;
            ++state;
            break;

        case(3):
            calculateMeans(ax, ay, az, gx, gy, gz, meanIterations);
            break;

        case(4): {
            gyroX->addData(t, gx_offset);
            gyroY->addData(t, gy_offset);
            gyroZ->addData(t, gz_offset);

            accX->addData(t, ax_offset);
            accY->addData(t, ay_offset);
            accZ->addData(t, az_offset);

            auto abs_ax = abs(mean_ax);
            auto abs_ay = abs(mean_ay);
            auto abs_az = abs(gravityConstant - mean_az);

            auto abs_gx = abs(mean_gx);
            auto abs_gy = abs(mean_gy);
            auto abs_gz = abs(mean_gz);

            qDebug() <<endl<<" gx_offset: "<< gx_offset <<" gy_offset: "<< gy_offset <<" gz_offset: "<< gz_offset <<" ax_offset: "<< ax_offset <<" ay_offset: "<< ay_offset <<" az_offset: "<< az_offset;
            qDebug() <<" gx_mean: "<< abs_gx <<" gy_mean: "<< abs_gy <<" gz_mean: "<< abs_gz <<" ax_mean: "<< abs_ax <<" ay_mean: "<< abs_ay <<" az_mean: "<< abs_az;

            auto sum = abs_ax + abs_ay + abs_az + abs_gx + abs_gy + abs_gz;
            if(sum < bestSum) {
                bestSum = sum;
                mean_ax_best = abs_ax; mean_ay_best = abs_ay; mean_az_best = abs_az; mean_gx_best = abs_gx; mean_gy_best = abs_gy; mean_gz_best = abs_gz;
                ax_offset_best = ax_offset; ay_offset_best = ay_offset; az_offset_best = az_offset; gx_offset_best = gx_offset; gy_offset_best = gy_offset; gz_offset_best = gz_offset;
            }


            int ready = 0;
            if (abs_ax <= acel_deadzone) ready++;
            else ax_offset = ax_offset - (mean_ax  * 0.5);
            if (abs_ay <= acel_deadzone) ready++;
            else ay_offset = ay_offset - (mean_ay  * 0.5);
            if (abs_az <= acel_deadzone) ready++;
            else az_offset = az_offset + ((gravityConstant - mean_az)  * 0.5);
            if (abs_gx <= gyro_deadzone) ready++;
            else gx_offset = gx_offset - (mean_gx  * 0.5);
            if (abs_gy <= gyro_deadzone) ready++;
            else gy_offset = gy_offset - (mean_gy  * 0.5);
            if (abs_gz <= gyro_deadzone) ready++;
            else gz_offset = gz_offset - (mean_gz  * 0.5);

            ++NOC;

            if (ready == 6 || NOC > 5000000)  {
                qDebug() <<endl<<" -------------------END ALGO---------------------";
                return;
            } else {
                buff_ax = 0.0f; buff_ay = 0.0f; buff_az = 0.0f; buff_gx = 0.0f; buff_gy = 0.0f; buff_gz = 0.0f;
                iterator = 0;
                state = 3;
            }
            break;
        }

        default:
            state = 0;
            break;
        }
    }
    qDebug() <<endl<<" -------------------END FILE---------------------";
    qDebug() <<"BEST PARAMETERS:";
    qDebug() <<endl<<" gx_offset: "<< gx_offset_best <<" gy_offset: "<< gy_offset_best <<" gz_offset: "<< gz_offset_best <<" ax_offset: "<< ax_offset_best <<" ay_offset: "<< ay_offset_best <<" az_offset: "<< az_offset_best;
    qDebug() <<" gx_mean: "<< mean_gx_best <<" gy_mean: "<< mean_gy_best <<" gz_mean: "<< mean_gz_best <<" ax_mean: "<< mean_ax_best <<" ay_mean: "<< mean_ay_best <<" az_mean: "<< mean_az_best;


    accX ->parentPlot()->replot();
    accY ->parentPlot()->replot();
    accZ ->parentPlot()->replot();

    gyroX ->parentPlot()->replot();
    gyroY ->parentPlot()->replot();
    gyroZ ->parentPlot()->replot();


    gyroFile.close();
    accFile.close();
}

void KalmanIMUTest::calculateMeans(double ax, double ay, double az, double gx, double gy, double gz, int meanIterations)
{
    if(iterator < (meanIterations + 101)) {
        if (iterator > 100 && iterator <= (meanIterations + 100)) { //First 100 measures are discarded
            buff_ax += ax;
            buff_ay += ay;
            buff_az += az;
            buff_gx += gx;
            buff_gy += gy;
            buff_gz += gz;
        }
        if (iterator == (meanIterations + 100)) {
            mean_ax = buff_ax / meanIterations;
            mean_ay = buff_ay / meanIterations;
            mean_az = buff_az / meanIterations;
            mean_gx = buff_gx / meanIterations;
            mean_gy = buff_gy / meanIterations;
            mean_gz = buff_gz / meanIterations;
        }
        ++iterator;
    } else {
        ++state;
    }
}
