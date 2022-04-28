#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    initPlot();
    kalmanTest = new KalmanIMUTest();

    // quaternion test kalman ---------------------------------------------------------------------------------------
    kalmanTest->testQuaternionKalman(alpha, beta, gamma, "Gyroscope.csv", "Accelerometer.csv");         // test kalman from files and plot to graph
    kalmanTest->csvXYZToGraph("Rotation.csv", alpha_real, beta_real, gamma_real, 10, (180.0f / M_PI));  // plot real angle data from dataset
    //---------------------------------------------------------------------------------------------------------------


    //kalmanTest->csvXYZToGraph("Gravity.csv", alpha, beta, gamma, 1, 1.0f);
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



