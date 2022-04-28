#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "kalmanImuTest.h"
#include <QMainWindow>
#include <qcustomplot.h>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void initPlot();

private:

    Ui::MainWindow *ui;

    QCPGraph* alpha;
    QCPGraph* beta;
    QCPGraph* gamma;

    QCPGraph* alpha_real;
    QCPGraph* beta_real;
    QCPGraph* gamma_real;

    KalmanIMUTest * kalmanTest;
};
#endif // MAINWINDOW_H
