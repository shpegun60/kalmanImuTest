#include "kalman_filter.h"
#include <stdlib.h>
#include "smart_assert/smart_assert.h"


KalmanFilter* kalmanCreate(initKalmanFunction userInit, unsigned int x, unsigned int z, unsigned int u)
{
    M_Assert_BreakSaveCheck((x == 0 || z == 0), "kalmanCreate: Give me positive values for dimensions genius", return NULL);

    KalmanFilter* m;
    m = (KalmanFilter*)malloc(sizeof(KalmanFilter));
    M_Assert_BreakSaveCheck((m == NULL), "kalmanCreate: no memory for allocation structure", return NULL);

    ////    * Result:
    m->X_est = matrixCreate(x, 1);  // estimated system state X[n,n]
    m->P_est = matrixCreate(x, x);  // estimated covariance P[n,n]

    ////    * Predict step:
    // state prediction
    m->X_pred = matrixCreate(x, 1);                             // X[n+1,n]
    m->F = matrixCreate(x, x);                                  // transition matrix F
    // drive
    if(u == 0) { // if not driving
        m->G = NULL;
        m->U = NULL;
        m->DriveMult = NULL;
    } else {
        m->G = matrixCreate(x, u);                                  // influence matrix G
        m->U = matrixCreate(u, 1);                                  // drive matrix U_n
        m->DriveMult = createResultMulMatrix(m->G, m->U);           // result multiplication ==> G * U_n in predict equation
    }

    // covariance prediction
    m->F_t = createResultTransMatrix(m->F);        // transition matrix transposed  F^T
    m->P_pred = matrixCreate(x, x);                // predict covariance matrix P[n+1,n]
    m->Q = matrixCreate(x, x);                     // system emulation covariations matrix Q_n

    ////    * Update step:

    // koefs
    m->K = matrixCreate(x, z);                     // Koefficients matrix K_n
    m->K_tmp = matrixCreate(x, z);                 // Matrix equal K_n to save multiplication (P[n,n-1] * H^T)
    m->S = matrixCreate(z, z);                     // matrix S_n
    m->S_inv = matrixCreate(z, z);                 // matrix S_n^-1
    m->R = matrixCreate(z, z);                     // measurments covariance matrix R_n
    // system state
    m->H = matrixCreate(z, x);                     // system measurments matrix H
    m->H_t = createResultTransMatrix(m->H);        // system measurments matrix H^T
    m->Z = matrixCreate(z, 1);                     // measurments matrix Z_n
    m->Update_derivative = matrixCreate(z, 1);     // multiplication result matrix ==> (Z_n − H * X[n,n−1]) in system update equaluation
    // covariance
    m->KHI = createResultMulMatrix(m->K, m->H);    // multiplication result matrix ==> (I − K_n * H)
    m->I = eye(x);                                 // identity matrix I for estimated covariance update

    if(userInit != NULL)  {
        userInit(m->X_est, m->P_est, m->F, m->F_t, m->G, m->Q, m->R, m->H, m->H_t);
    }

    return m;
}


int kalmanPredict(KalmanFilter* m)
{
    M_Assert_Break((m == NULL), "kalmanPredict: m is not exists", return KALMAN_ERR);
    /*
    *****************************************
    * Predict step:
    * 1) X[n+1,n] = F * X[n,n] + G * U_n;       // predict system
    * 2) P[n+1,n] = F * P[n,n] * F^T + Q_n;     // predict covariance
    *****************************************
    */

    // 1)
    multiply(m->F, m->X_est, m->X_pred);        // F * X[n,n] = X[n+1,n]                // 200
    multiply(m->G, m->U, m->DriveMult);         // G * U_n = DriveMult                  // 60
    add(m->X_pred, m->DriveMult, m->X_pred);    // X[n+1,n] = X[n+1,n] + DriveMult      // 10

    // 2)
    multiply(m->P_est, m->F_t, m->KHI);         // P[n,n] * F^T = KHI                   //2000
    multiply(m->F, m->KHI, m->P_pred);          // F * KHI = P[n+1,n]                   //2000
    add(m->P_pred, m->Q, m->P_pred);            // P[n+1,n] = P[n+1,n] + Q_n            //100
    // total floating operation: 4370

    return KALMAN_OK;
}

int kalmanPredict_withoutDrive(KalmanFilter* m)
{
    M_Assert_Break((m == NULL), "kalmanPredict: m is not exists", return KALMAN_ERR);
    /*
    *****************************************
    * Predict step:
    * 1) X[n+1,n] = F * X[n,n]                  // predict system
    * 2) P[n+1,n] = F * P[n,n] * F^T + Q_n;     // predict covariance
    *****************************************
    */

    // 1)
    multiply(m->F, m->X_est, m->X_pred);        // F * X[n,n] = X[n+1,n]

    // 2)
    multiply(m->P_est, m->F_t, m->KHI);         // P[n,n] * F^T = KHI
    multiply(m->F, m->KHI, m->P_pred);          // F * KHI = P[n+1,n]
    add(m->P_pred, m->Q, m->P_pred);            // P[n+1,n] = P[n+1,n] + Q_n

    return KALMAN_OK;
}

int kalmanUpdate(KalmanFilter* m)
{
    M_Assert_Break((m == NULL), "kalmanUpdate: m is not exists", return KALMAN_ERR);
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

    // 3)
    multiply(m->P_pred, m->H_t, m->K_tmp);          // (P[n,n-1] * H^T) = K_tmp                                     //800
    multiply(m->H, m->K_tmp, m->S);                 // H * K_tmp = S_n                                              // 320
    add(m->S, m->R, m->S);                          // S_n = S_n + R_n                                              // 16

    // 4)
    gluInvertMatrix4x4_fastest(m->S, m->S_inv);     // S_n ==> S_n^-1                                               // 200
    multiply(m->K_tmp, m->S_inv, m->K);             // K_tmp * S_n^-1 = K_n                                         // 320

    // 5)
    multiply(m->H, m->X_pred, m->Update_derivative);            // H * X[n,n−1] = Update_derivative                 // 80
    sub(m->Z, m->Update_derivative, m->Update_derivative);      // Z_n − Update_derivative = Update_derivative      // 4
    multiply(m->K, m->Update_derivative, m->DriveMult);         // K_n * Update_derivative = DriveMult              // 80
    add(m->X_pred, m->DriveMult, m->X_est);                     // X[n,n] = X[n,n−1] + DriveMult

    // 6)
    multiply(m->K, m->H, m->KHI);          // K_n * H = KHI                                                         // 800
    sub(m->I, m->KHI, m->KHI);             // (I − KHI) = KHI                                                       // 100
    multiply(m->KHI, m->P_pred, m->P_est); // KHI * P[n,n−1] = P[n,n]                                               // 2000
    // total floating operation: 4720
    return KALMAN_OK;
}
// total operation Kalman: 9090 ~ 10000

void printkalman(KalmanFilter* m)
{
    M_Assert_BreakSaveCheck((m == NULL), "printkalman: kalman is not exist", return);
    showmat(m->X_est, (char *)"X_est:");
    showmat(m->P_est, (char *)"P_est:");
    showmat(m->X_pred, (char *)"X_pred:");
    showmat(m->F, (char *)"F:");
    showmat(m->G, (char *)"G:");
    showmat(m->U, (char *)"U:");
    showmat(m->DriveMult, (char *)"DriveMult:");
    showmat(m->F_t, (char *)"F_t:");
    showmat(m->P_pred, (char *)"P_pred:");
    showmat(m->Q, (char *)"Q:");
    showmat(m->K, (char *)"K:");
    showmat(m->K_tmp, (char *)"K_tmp:");
    showmat(m->S, (char *)"S:");
    showmat(m->S_inv, (char *)"S_inv:");
    showmat(m->R, (char *)"R:");
    showmat(m->H, (char *)"H:");
    showmat(m->H_t, (char *)"H_t:");
    showmat(m->Z, (char *)"Z:");
    showmat(m->Update_derivative, (char *)"Update_derivative:");
    showmat(m->KHI, (char *)"KHI:");
    showmat(m->I, (char *)"I:");
}
