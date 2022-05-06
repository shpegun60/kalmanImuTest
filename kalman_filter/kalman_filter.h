/**
 * @file    KALMAN_FILTER.h
 * @brief
 * @date
 */

#ifndef __KALMAN_FILTER
#define __KALMAN_FILTER

#include "kalman_port.h"
#include "matrix.h"

/**
 * Data structure to hold a kalman filter system state.
 */

typedef void (*initKalmanFunction)(
        Mat* X_est, // init X[n,n]
        Mat* P_est, // init P[n,n]
        Mat* F,     // init transition state matrix
        Mat* F_t,   // init transition state matrix transposed
        Mat* G,     // init drive matrix
        Mat* Q,     // init emulation covariance matrix
        Mat* R,     // init measurment covariance matrix
        Mat* H,     // init Observation Matrix
        Mat* H_t    // init Observation Matrix transposed
        );

typedef void (*updateKalmanInputData) (
        Mat* Z,     // update input measurement matrix
        Mat* F,     // update transition state matrix
        Mat* F_t,   // update transition state matrix transposed
        Mat* G,     // update drive matrix
        Mat* U,     // update input drive measurments
        Mat* Q,     // update emulation covariance matrix
        Mat* R      // update measurment covariance matrix
        );

typedef struct {
    /*
    *****************************************
    * Result
    *****************************************
    */

    Mat* X_est; // estimated system state X[n,n]
    Mat* P_est; // estimated covariance P[n,n]

    /*
    *****************************************
    * Predict step:
    * 1) X[n+1,n] = F * X[n,n] + G * U_n;   // predict system
    * 2) P[n+1,n] = F * P[n,n] * F^T + Q_n; // predict covariance
    *****************************************
    */

    // state prediction
    Mat* X_pred;                // X[n+1,n]
    Mat* F;                     // transition matrix F // USER OWERWRITE
    // drive
    Mat* G;                     // influence matrix G // USER OWERWRITE
    Mat* U;                     // drive matrix U_n   // USER OWERWRITE
    Mat* DriveMultPredict;      // result multiplication ==> G * U_n in predict equation
    // covariance prediction
    Mat* F_t;                   // transition matrix transposed  F^T // USER OWERWRITE
    Mat* P_pred;                // predict covariance matrix P[n+1,n]
    Mat* Q;                     // system emulation covariations matrix Q_n // USER OWERWRITE

    /*
    *****************************************
    * Update step:
    * 3) S_n = H * P[n,n-1] * H^T + R_n;                   // update S_n matrix
    * 4) K_n = (P[n,n-1] * H^T) / S_n;                     // update K_n matrix, based on S_n
    * 5) X[n,n] = X[n,n−1] + K_n * (Z_n − H * X[n,n−1]);   // update result system X[n,n]
    *
    *          - > P[n,n] = (I − K_n * H) * P[n,n−1];      // simple estimated covariance
    * 6) (or)  |
    *          - > P[n,n] = (I − K_n * H) * P[n,n−1] * (I − K_n * H)^T + K_n * R_n * K_n^T; // extended covariance
    * 7) go to 1);
    *
    *****************************************
    */

    // koefs
    Mat* K;                     // Koefficients matrix K_n
    Mat* K_tmp;                 // Matrix equal K_n to save multiplication (P[n,n-1] * H^T)
    Mat* S;                     // matrix S_n
    Mat* S_inv;                 // matrix S_n^-1
    Mat* R;                     // measurments covariance matrix R_n // USER OWERWRITE
    // system state
    Mat* H;                     // system measurments matrix H
    Mat* H_t;                   // system measurments matrix H^T
    Mat* Z;                     // measurments matrix Z_n // USER OWERWRITE
    Mat* Update_derivative;     // multiplication result matrix ==> (Z_n − H * X[n,n−1]) in system update equaluation
    Mat* DriveMultUpdate;       //buffer to equaluation update system K_n * (Z_n − H * X[n,n−1])
    // covariance
    Mat* KHI;                   // multiplication result matrix ==> (I − K_n * H)
    Mat* I;                     // identity matrix I for estimated covariance update

} KalmanFilter;


KalmanFilter* kalmanCreate(initKalmanFunction userInit, unsigned int x, unsigned int z, unsigned int u);
int kalmanPredict(KalmanFilter* m);
int kalmanPredict_withoutDrive(KalmanFilter* m);
int kalmanUpdate(KalmanFilter* m);

void printkalman(KalmanFilter* m);

#endif /* __KALMAN_FILTER */

