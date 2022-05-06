#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    initPlot();
    kalmanTest = new KalmanIMUTest();

//    //        // quaternion test kalman ---------------------------------------------------------------------------------------
            kalmanTest->testQuaternionKalman(alpha, beta, gamma, "Gyroscope.csv", "Accelerometer.csv", ANGLE_DATA);         // test kalman from files and plot to graph
            kalmanTest->csvXYZToGraph("Rotation.csv", alpha_real, beta_real, gamma_real, 10, (180.0f / M_PI));  // plot real angle data from dataset
    //        //---------------------------------------------------------------------------------------------------------------


//    //    // gravity test ---------------------------------------------------------------------------------------
//    kalmanTest->testQuaternionKalman(alpha, beta, gamma, "Gyroscope.csv", "Accelerometer.csv", GRAVITY_DATA);
//    kalmanTest->csvXYZToGraph("Accelerometer.csv", alpha_real, beta_real, gamma_real, 1, 1.0f);
    //    //---------------------------------------------------------------------------------------------------------------


//    //    // linear acceleration test ---------------------------------------------------------------------------------------
//    kalmanTest->testQuaternionKalman(alpha, beta, gamma, "Gyroscope.csv", "Accelerometer.csv", LINEAR_ACCERARATION_DATA);
//    kalmanTest->csvXYZToGraph("Magnetometer.csv", alpha_real, beta_real, gamma_real, 1, 1.0f);
    //    //---------------------------------------------------------------------------------------------------------------

    /// calibration accel & gyro from file
    //kalmanTest->calibrateAccelGyroFromFile(alpha, beta, gamma, alpha_real, beta_real, gamma_real, "Gyroscope.csv", "Accelerometer.csv", 9.8f, 200, 0.00001, 0.00001);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initPlot()
{
    // add test graph
    alpha = ui->qplot1->addGraph();
    beta = ui->qplot1->addGraph();
    gamma = ui->qplot1->addGraph();

    alpha->setPen(QPen(Qt::blue));
    beta->setPen(QPen(Qt::red));
    gamma->setPen(QPen(Qt::green));

    alpha->setScatterStyle(QCPScatterStyle::ssCircle);
    beta->setScatterStyle(QCPScatterStyle::ssCircle);
    gamma->setScatterStyle(QCPScatterStyle::ssCircle);

    ui->qplot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->qplot1->legend->setVisible(true);

    // add real graph
    alpha_real = ui->qplot2->addGraph();
    beta_real = ui->qplot2->addGraph();
    gamma_real = ui->qplot2->addGraph();

    alpha_real->setPen(QPen(Qt::blue));
    beta_real->setPen(QPen(Qt::red));
    gamma_real->setPen(QPen(Qt::green));


    alpha_real->setScatterStyle(QCPScatterStyle::ssCircle);
    beta_real->setScatterStyle(QCPScatterStyle::ssCircle);
    gamma_real->setScatterStyle(QCPScatterStyle::ssCircle);

    ui->qplot2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->qplot2->legend->setVisible(true);
}



