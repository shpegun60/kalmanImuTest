#include "imu.h"
#include <stdlib.h>
#include "smart_assert/smart_assert.h"
#include <assert.h>


#define INV_SQRT_2 0.70710678118654752440f
#define EPSILON_QUATER 0.9f

#define FABS(x) ((x) < 0.0f ? -(x) : (x))
#define TO_RAD(x) (((x) / 180.0f) * M_PI)
#define TO_DEG(x) (((x) * 180.0f) / M_PI)
static_assert (__builtin_types_compatible_p(MAT_TYPE, float), "IMU: MATH_TYPE must be float or rewrite /void Quaternion_multiply_to_arrayLN(Quaternion* q1, Quaternion* q2, float** output)/ in quaternion lib to /void Quaternion_multiply_to_arrayLN(Quaternion* q1, Quaternion* q2, USER_TYPE** output)/ and rewrite this assert to USER_TYPE");

/*
 *
 *           | x | xh| y | yh| z | zh| q0| q1| q1| q3|
 *        ---|---|---|---|---|---|---|---|---|---|---|
 * F = [   x | 1 | dt| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0
 *        xh | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1
 *         y | 0 | 0 | 1 | dt| 0 | 0 | 0 | 0 | 0 | 0 | 2
 *        yh | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 3
 *        z  | 0 | 0 | 0 | 0 | 1 | dt| 0 | 0 | 0 | 0 | 4
 *        zh | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 5
 *        q0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |-Tx|-Ty|-Tz| 6
 *        q1 | 0 | 0 | 0 | 0 | 0 | 0 | Tx| 1 | Tz|-Ty| 7
 *        q2 | 0 | 0 | 0 | 0 | 0 | 0 | Ty|-Tz| 1 | Tx| 8
 *        q3 | 0 | 0 | 0 | 0 | 0 | 0 | Tz| Ty|-Tx| 1 | 9    ]
 *----------------------------------------------------
 *             0   1   2   3   4   5   6   7   8   9
 *
 *  Tx = ( dt/2 ) * w_x
 *  Ty = ( dt/2 ) * w_y
 *  Tz = ( dt/2 ) * w_z
 *
 *
 *           |xhh|yhh|zhh|
 *        ---|---|---|---|
 * G = [   x |Tg | 0 | 0 | 0
 *        xh |dt | 0 | 0 | 1
 *         y | 0 |Tg | 0 | 2
 *        yh | 0 |dt | 0 | 3
 *        z  | 0 | 0 |Tg | 4
 *        zh | 0 | 0 |dt | 5
 *        q0 | 0 | 0 | 0 | 6
 *        q1 | 0 | 0 | 0 | 7
 *        q2 | 0 | 0 | 0 | 8
 *        q3 | 0 | 0 | 0 | 9    ]
 *------------------------
 *             0   1   2
 *
 * Tg = ( dt^2 ) / 2
 *
 *
 * H = [    | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
 *          | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
 *          | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
 *          | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |      ]
 *
 */


//typedef struct {
//    // variances
//    float gyroVarianceVect[3];   // gyroscope variance vector [x, y, z]
//    float accVarianceVect[3];    // accelerometer variance vector [x, y, z]
//    float magVarianceVect[3];    // magnetometer variance vector [x, y, z]

//} IMUInit_struct;

#define KALx 10  // x_k system state
#define KALz 4   // z_k measurment
#define KALu 3   // u_k drive

IMU10Dof* imuCreate(IMUInit_struct* init)
{
    M_Assert_BreakSaveCheck((init == NULL), "imuCreate:INIT is nullptr", return NULL);
    M_Assert_BreakSaveCheck((init->accConst_u < 0.0f || init->accConst_u > 1.0f), "imuCreate: value of u out of range", return NULL);
    IMU10Dof* imu = (IMU10Dof*)malloc(sizeof(IMU10Dof));
    M_Assert_BreakSaveCheck((imu == NULL), "imuCreate: no memory for allocation structure", return NULL);

    imu->kalman = kalmanCreate(NULL, KALx, KALz, KALu);

    /*
     * ******************************************************************************************
     * KALMAN MATRIX user init functions
     *******************************************************************************************
     *
    */

    // F, F_t init ------------------------------------------------------------------------
    float dt = init->dt_init_sec;
    float F_init [KALx][KALx] =
    {
        {1, dt, 0,  0,  0,  0,  0,  0,  0,  0},
        {0, 1,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  1,  dt, 0,  0,  0,  0,  0,  0},
        {0, 0,  0,  1,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  1,  dt, 0,  0,  0,  0},
        {0, 0,  0,  0,  0,  1,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  1,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  1,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  1,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  1}
    };
    matrixInitFromArr(imu->kalman->F, F_init[0]);
    matrixInitFromArr_T(imu->kalman->F_t, F_init[0]);

    // X_est init ------------------------------------------------------------------------
    float X_est_init[KALx][1] =
    {
        {0},
        {0},
        {0},
        {0},
        {0},
        {0},
        {1},
        {0},
        {0},
        {0}
    };
    matrixInitFromArr(imu->kalman->X_est, X_est_init[0]);

    // P_est init ------------------------------------------------------------------------
    float P_init[KALx][KALx] =
    {
        {500, 0,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 500,  0,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  500,  0,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  500,  0,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  500,  0,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  500,  0,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  500,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  500,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  500,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  500}
    };
    matrixInitFromArr(imu->kalman->P_est, P_init[0]);

    // G init ------------------------------------------------------------------------
    float Tg = ( dt * dt ) / 2.0f;
    float G_init [KALx][KALu] =
    {
        {Tg, 0,  0   },
        {dt, 0,  0   },
        {0,  Tg, 0   },
        {0,  dt, 0,  },
        {0,  0,  Tg, },
        {0,  0,  dt, },
        {0,  0,  0,  },
        {0,  0,  0,  },
        {0,  0,  0,  },
        {0,  0,  0,  }
    };
    matrixInitFromArr(imu->kalman->G, G_init[0]);

    // H init ------------------------------------------------------------------------
    float H_init[KALz][KALx] =
    {
        {0, 0,  0,  0,  0,  0,  1,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  1,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  1,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  1}
    };
    matrixInitFromArr(imu->kalman->H, H_init[0]);
    matrixInitFromArr_T(imu->kalman->H_t, H_init[0]);


    // Q init ------------------------------------------------------------------------
    float sQR = 1e-6;
    float sQx = init->accVarianceVect[0];
    float sQy = init->accVarianceVect[1];
    float sQz = init->accVarianceVect[2];

    float T11 = (dt * dt * dt * dt) / 4.0;
    float T12 = (dt * dt * dt) / 2.0;
    float T21 = (dt * dt * dt) / 2.0;
    float T22 = (dt * dt);

    float Q_init[KALx][KALx] =
    {
        {T11 * sQx,     T12 * sQx,              0,                      0,                              0,                          0,                  0,  0,  0,  0},
        {T21 * sQx,     T22 * sQx,              0,                      0,                              0,                          0,                  0,  0,  0,  0},
        {0,             0,                      T11 * sQy,              T12 * sQy,                      0,                          0,                  0,  0,  0,  0},
        {0,             0,                      T21 * sQy,              T22 * sQy,                      0,                          0,                  0,  0,  0,  0},
        {0,             0,                      0,                      0,                              T11 * sQz,                  T12 * sQz,          0,  0,  0,  0},
        {0,             0,                      0,                      0,                              T21 * sQz,                  T22 * sQz,          0,  0,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  sQR,  0,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  sQR,  0,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  0,  sQR,  0},
        {0,             0,                      0,                      0,                              0,                          0,                  0,  0,  0,  sQR}
    };
    matrixInitFromArr(imu->kalman->Q, Q_init[0]);

    // R init ------------------------------------------------------------------------
    float sR = 0.0015f;
    float R_init[KALz][KALz] =
    {
        {sR, 0,  0,  0},
        {0, sR,  0,  0},
        {0, 0,  sR,  0},
        {0, 0,  0,  sR}
    };
    matrixInitFromArr(imu->kalman->R, R_init[0]);


    /*
     * ******************************************************************************************
     * IMU user init functions
     *******************************************************************************************
     *
    */

    if(init->coordinateType == NED) {
        imu->accDeltaQuaterFinder = accDeltaQuaterFinder_NED;
        imu->magDeltaQuaterFinder = magDeltaQuaterFinder_NED;
    } else {
        imu->accDeltaQuaterFinder = accDeltaQuaterFinder_ENU;
        imu->magDeltaQuaterFinder = magDeltaQuaterFinder_ENU;
    }

    imu->accConst_u = init->accConst_u;
    imu->magConst_u = init->magConst_u;
    for(unsigned i = 0; i < 3; ++i) {
        imu->grav[i] = 0.0f;
        imu->rotateAxisErr[i] = 0.0f;
        imu->grav_norm[i] = 0.0f;

        if(init->coordinateType == NED) {
            imu->calibrationGrav[i] = -init->gravityConstVect[i];
        } else {
            imu->calibrationGrav[i] = init->gravityConstVect[i];
        }

        imu->accelBiasVect[i] = init->accelBiasVect[i];
        imu->gyroBiasVect[i] = init->gyroBiasVect[i];
    }

    Quaternion_setIdentity(&imu->q_ae);
    Quaternion_setIdentity(&imu->q_a);
    Quaternion_setIdentity(&imu->q_i);

    imu->Tx = 0.0f;
    imu->Ty = 0.0f;
    imu->Tz = 0.0f;
    imu->halfT = 0.0f;
    imu->recipNorm = 0.0f;


    // for update Q ----------------------------------------------
    imu->Q_quat = matrixCreate(4, 4);
    imu->G = matrixCreate(4, 3);
    imu->G_t = matrixCreate(3, 4);
    imu->Noise = matrixCreate(3, 3);
    imu->NOISE_RES = createResultMulMatrix(imu->Noise, imu->G_t);

    for(unsigned i = 0; i < imu->Noise->col; ++i) {
        imu->Noise->data[i][i] = init->gyroVarianceVect[i];
    }

    imu->constT_4 = 0.0f;

    // R
    imu->J = matrixCreate(4, 3);
    imu->J_t = matrixCreate(3, 4);
    imu->Noise_acc = matrixCreate(3, 3);
    imu->NOISE_R_RES = createResultMulMatrix(imu->Noise_acc, imu->J_t);

    return imu;
}


__attribute__((unused)) static void computeGravVector(const float q0, const float q1, const float q2, const float q3, float* const grav, const float* const constant)
{
    M_Assert_Break((grav == NULL || constant == NULL), "computeGravVector: grav vector is null pointer", return);
    grav[0] = 2.0f * ((q1 * q3) - (q0 * q2)) * constant[0];
    grav[1] = 2.0f * ((q0 * q1) + (q2 * q3)) * constant[1];
    grav[2] = ((q0 * q0) - (q1 * q1) - (q2 *q2) + (q3 * q3)) * constant[2];
}

__attribute__((unused)) static void vecMulCrossAngle(const float* const a, const float* const b, float* const r, float* const angle)
{
    M_Assert_Break((a == NULL || b == NULL || r == NULL || angle == NULL), "vecMulCross: vectors is null ptr", return);

    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];

    // very precision method with atan2
    float dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    float det = fastSqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    *angle = atan2(det, dot);
}

__attribute__((unused)) static float accelCalcPitch(float x, float y, float z)
{
    return atan2(y, fastSqrt(x * x + z * z));
}

__attribute__((unused)) static float accelCalcRoll(float x, float y, float z)
{
    return -atan2(x, fastSqrt(y * y + z * z));
}

int imuProceed(IMU10Dof* imu, IMUinput * data) // original kalman fusion
{
    M_Assert_Break((imu == NULL || data == NULL), "imuProceed: vectors is null ptr", return KALMAN_ERR);
    M_Assert_Break((imu->accDeltaQuaterFinder == NULL || imu->magDeltaQuaterFinder == NULL), "imuProceed: NULL delta quater functions", return KALMAN_ERR);

    /*
     * ***************************************************
     * Predict step 1
     * ***************************************************
     */

    // update U
    for (unsigned int i  = 0; i < 3; ++i) {
        data->a[i] += imu->accelBiasVect[i];
        data->g[i] += imu->gyroBiasVect[i];
        imu->kalman->U->data[i][0] = data->a[i] - imu->grav[i];
        //imu->kalman->U->data[i][0] = (FABS(imu->kalman->U->data[i][0]) < 0.001f) ? 0.0f : imu->kalman->U->data[i][0];
    }

    // update F, F_t
    imu->kalman->F_t->data[5][4] = imu->kalman->F->data[4][5] = imu->kalman->F_t->data[3][2]
            = imu->kalman->F->data[2][3] = imu->kalman->F_t->data[1][0] = imu->kalman->F->data[0][1]
            = data->dt_sec;

    imu->halfT = 0.5f * data->dt_sec;
    imu->Tx = imu->halfT * data->g[0];
    imu->Ty = imu->halfT * data->g[1];
    imu->Tz = imu->halfT * data->g[2];

    imu->kalman->F_t->data[6][7] = imu->kalman->F->data[7][6] = imu->kalman->F_t->data[9][8] = imu->kalman->F->data[8][9] = imu->Tx;
    imu->kalman->F_t->data[6][8] = imu->kalman->F->data[8][6] = imu->kalman->F_t->data[7][9] = imu->kalman->F->data[9][7] = imu->Ty;
    imu->kalman->F_t->data[8][7] = imu->kalman->F->data[7][8] = imu->kalman->F_t->data[6][9] = imu->kalman->F->data[9][6] = imu->Tz;

    imu->kalman->F_t->data[7][6] = imu->kalman->F->data[6][7] = imu->kalman->F_t->data[8][9] = imu->kalman->F->data[9][8] = -imu->Tx;
    imu->kalman->F_t->data[8][6] = imu->kalman->F->data[6][8] = imu->kalman->F_t->data[9][7] = imu->kalman->F->data[7][9] = -imu->Ty;
    imu->kalman->F_t->data[9][6] = imu->kalman->F->data[6][9] = imu->kalman->F_t->data[7][8] = imu->kalman->F->data[8][7] = -imu->Tz;


    // update G
    imu->kalman->G->data[4][2] = imu->kalman->G->data[2][1] = imu->kalman->G->data[0][0] = (data->dt_sec * data->dt_sec) * 0.5f;
    imu->kalman->G->data[5][2] = imu->kalman->G->data[3][1] = imu->kalman->G->data[1][0] = data->dt_sec;


    //Update Q
    imu->G_t->data[2][2] = imu->G->data[2][2] = imu->G_t->data[0][0] = imu->G->data[0][0] = imu->kalman->X_est->data[7][0]; // q1
    imu->G_t->data[0][3] = imu->G->data[3][0] = imu->G_t->data[1][0] = imu->G->data[0][1] = imu->kalman->X_est->data[8][0]; // q2
    imu->G_t->data[1][1] = imu->G->data[1][1] = imu->G_t->data[2][0] = imu->G->data[0][2] = imu->kalman->X_est->data[9][0]; // q3

    imu->G_t->data[2][3] = imu->G->data[3][2] = imu->G_t->data[1][2] = imu->G->data[2][1] = imu->G_t->data[0][1] = imu->G->data[1][0] = -imu->kalman->X_est->data[6][0]; // -q0
    imu->G_t->data[1][3] = imu->G->data[3][1] = -imu->kalman->X_est->data[7][0]; // -q1
    imu->G_t->data[2][1] = imu->G->data[1][2] = -imu->kalman->X_est->data[8][0]; // -q2
    imu->G_t->data[0][2] = imu->G->data[2][0] = -imu->kalman->X_est->data[9][0]; // -q3
    imu->constT_4 = (data->dt_sec * data->dt_sec) * 0.25f;

    multiply(imu->Noise, imu->G_t, imu->NOISE_RES);
    multiply(imu->G, imu->NOISE_RES, imu->Q_quat);
    scalarmultiply(imu->Q_quat, imu->Q_quat, imu->constT_4);

    for(unsigned int i = 0; i < imu->Q_quat->row; ++i) {
        for(unsigned int j = 0; j < imu->Q_quat->col; ++j) {
            imu->kalman->Q->data[6 + i][6 + j] = imu->Q_quat->data[i][j];
        }
    }

    kalmanPredict(imu->kalman);
    //matrixCopy(imu->kalman->X_pred, imu->kalman->X_est);
    Quaternion_normalize_vect(&imu->kalman->X_pred->data[6][0], &imu->kalman->X_pred->data[7][0], &imu->kalman->X_pred->data[8][0], &imu->kalman->X_pred->data[9][0]);

    /*
     * ***************************************************
     * Update Step 2
     * ***************************************************
     */

    ////////    //    // normalize accel
    imu->recipNorm = invSqrt(data->a[0] * data->a[0] + data->a[1] * data->a[1] + data->a[2] * data->a[2]);
    data->a[0] *= imu->recipNorm;
    data->a[1] *= imu->recipNorm;
    data->a[2] *= imu->recipNorm;


    Quaternion_set(imu->kalman->X_pred->data[6][0], imu->kalman->X_pred->data[7][0], imu->kalman->X_pred->data[8][0], imu->kalman->X_pred->data[9][0], &imu->q_a); // qk
    Quaternion_rotate(&imu->q_a, data->a, imu->rotateAxisErr);
    imu->accDeltaQuaterFinder(imu, data);

    // filtration
    if(imu->q_ae.w > EPSILON_QUATER) {
        // LERP
        Quaternion_lerp(&imu->q_i, &imu->q_ae, imu->accConst_u, &imu->q_ae);
    } else {
        //SLERP
        Quaternion_slerp(&imu->q_i, &imu->q_ae, imu->accConst_u, &imu->q_ae);
    }

    Quaternion_multiply(&imu->q_a, &imu->q_ae, &imu->q_a);
    //Quaternion_multiply_to_arrayLN(&imu->q_a, &imu->q_ae, imu->kalman->Z->data);
    //Quaternion_multiply(&imu->q_a, &imu->q_ae, &imu->RES);

    ////////    //    // normalize magnetometer
    imu->recipNorm = invSqrt(data->m[0] * data->m[0] + data->m[1] * data->m[1] + data->m[2] * data->m[2]);
    data->m[0] *= imu->recipNorm;
    data->m[1] *= imu->recipNorm;
    data->m[2] *= imu->recipNorm;


    Quaternion_rotate(&imu->q_a, data->m, imu->rotateAxisErr);

    imu->magDeltaQuaterFinder(imu, data);

    // filtration
    if(imu->q_ae.w > EPSILON_QUATER) {
        // LERP
        Quaternion_lerp(&imu->q_i, &imu->q_ae, imu->magConst_u, &imu->q_ae);
    } else {
        //SLERP
        Quaternion_slerp(&imu->q_i, &imu->q_ae, imu->magConst_u, &imu->q_ae);
    }

    Quaternion_multiply_to_arrayLN(&imu->q_a, &imu->q_ae, imu->kalman->Z->data);
    //Quaternion_multiply(&imu->q_a, &imu->q_ae, &imu->RES);


    kalmanUpdate(imu->kalman);

    Quaternion_normalize_vect(&imu->kalman->X_est->data[6][0], &imu->kalman->X_est->data[7][0], &imu->kalman->X_est->data[8][0], &imu->kalman->X_est->data[9][0]);
    computeGravVector(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], imu->grav, imu->calibrationGrav); // a_calculated
    //computeGravVector(imu->kalman->Z->data[0][0], imu->kalman->Z->data[1][0], imu->kalman->Z->data[2][0], imu->kalman->Z->data[3][0], imu->grav, imu->calibrationGrav); // a_calculated
    Quaternion_set(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], &imu->q_a); // qk
    return KALMAN_OK;
}

// accelerometer delta quaternion finder for all coordinates --------------------------------------------
void accDeltaQuaterFinder_NED(IMU10Dof* imu, IMUinput * data) // NED gravity vector G = [0, 0, -1]
{
    //-g system
    imu->recipNorm = invSqrt(2.0f * (1.0f - imu->rotateAxisErr[2]));
    imu->q_ae.w    = (INV_SQRT_2 * fastSqrt(1.0f - imu->rotateAxisErr[2]));
    imu->q_ae.v[0] = -imu->rotateAxisErr[1] * imu->recipNorm;
    imu->q_ae.v[1] = imu->rotateAxisErr[0] * imu->recipNorm;
    imu->q_ae.v[2] = 0.0f;
    (void)data;
}

void accDeltaQuaterFinder_ENU(IMU10Dof* imu, IMUinput * data) // ENU gravity vector G = [0, 0,  1]
{
    //g system
    imu->recipNorm = invSqrt(2.0f * (imu->rotateAxisErr[2] + 1.0f));
    imu->q_ae.w    = (INV_SQRT_2 * fastSqrt(imu->rotateAxisErr[2] + 1.0f));
    imu->q_ae.v[0] = imu->rotateAxisErr[1] * imu->recipNorm;
    imu->q_ae.v[1] = -imu->rotateAxisErr[0] * imu->recipNorm;
    imu->q_ae.v[2] = 0.0f;
    (void)data;
}

// magnetometer delta quaternion finder for all coordinates --------------------------------------------
void magDeltaQuaterFinder_NED(IMU10Dof* imu, IMUinput * data) // NED geomagnetic field vector R = [mN, 0, mD]
{
    float G = imu->rotateAxisErr[0]*imu->rotateAxisErr[0] + imu->rotateAxisErr[1]*imu->rotateAxisErr[1];
    float G_k = fastSqrt(G);
    float tmp = fastSqrt(G + (imu->rotateAxisErr[0] * G_k));

    imu->q_ae.w    = tmp * invSqrt(2.0f * G);
    imu->q_ae.v[0] = 0.0f;
    imu->q_ae.v[1] = 0.0f;
    imu->q_ae.v[2] = -imu->rotateAxisErr[1] * INV_SQRT_2 * invSqrt(tmp * tmp);

    //    // hard version
    //        float G = imu->rotateAxisErr[0]*imu->rotateAxisErr[0] + imu->rotateAxisErr[1]*imu->rotateAxisErr[1] + 0.0001;
    //        imu->q_ae.w    = -imu->rotateAxisErr[1] * pow(G, -0.1e1 / 0.2e1) * pow(-0.2e1 * (sqrt(G) * imu->rotateAxisErr[0] - G) / G, -0.1e1 / 0.2e1);
    //        imu->q_ae.v[0] = 0.0f;
    //        imu->q_ae.v[1] = 0.0f;
    //        imu->q_ae.v[2] = sqrt(-0.2e1 * (sqrt(G) * imu->rotateAxisErr[0] - G) / G) / 0.2e1;

    (void)data;
}

void magDeltaQuaterFinder_ENU(IMU10Dof* imu, IMUinput * data) // ENU geomagnetic field vector R = [0, mN, -mD]
{
    float G = imu->rotateAxisErr[0]*imu->rotateAxisErr[0] + imu->rotateAxisErr[1]*imu->rotateAxisErr[1];
    float G_k = fastSqrt(G);
    float G_k_inv = invSqrt(G);

    imu->q_ae.v[2] = fastSqrt(G - imu->rotateAxisErr[1] * G_k) * invSqrt(2.0f * G);
    imu->q_ae.w    = imu->rotateAxisErr[0] * invSqrt(imu->q_ae.v[2] * imu->q_ae.v[2]) * 0.5f * G_k_inv;
    imu->q_ae.v[0] = 0.0f;
    imu->q_ae.v[1] = 0.0f;

    //    // hard version
    //    float G = imu->rotateAxisErr[0]*imu->rotateAxisErr[0] + imu->rotateAxisErr[1]*imu->rotateAxisErr[1];
    //    imu->q_ae.v[2] = sqrt(-0.2e1 * (imu->rotateAxisErr[1] * sqrt(G) - G) / G) / 0.2e1;
    //    imu->q_ae.w    = imu->rotateAxisErr[0] * pow(G, -0.1e1 / 0.2e1) * pow(-0.2e1 * (imu->rotateAxisErr[1] * sqrt(G) - G) / G, -0.1e1 / 0.2e1);
    //    imu->q_ae.v[0] = 0.0f;
    //    imu->q_ae.v[1] = 0.0f;

    (void)data;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//int imuProceed(IMU10Dof* imu, IMUinput * data) // original kalman fusion
//{
//    M_Assert_Break((imu == NULL || data == NULL), "imuProceed: vectors is null ptr", return KALMAN_ERR);
//    /*
//     * ***************************************************
//     * Predict step 1
//     * ***************************************************
//     */

//    // update U
//    for (unsigned int i  = 0; i < 3; ++i) {
//        data->a[i] += imu->accelBiasVect[i];
//        data->g[i] += imu->gyroBiasVect[i];
//        imu->kalman->U->data[i][0] = data->a[i] - imu->grav[i];
//        //imu->kalman->U->data[i][0] = (FABS(imu->kalman->U->data[i][0]) < 0.001f) ? 0.0f : imu->kalman->U->data[i][0];
//    }

//    // update F, F_t
//    imu->kalman->F_t->data[5][4] = imu->kalman->F->data[4][5] = imu->kalman->F_t->data[3][2]
//            = imu->kalman->F->data[2][3] = imu->kalman->F_t->data[1][0] = imu->kalman->F->data[0][1]
//            = data->dt;

//    imu->halfT = 0.5f * data->dt;
//    imu->Tx = imu->halfT * data->g[0];
//    imu->Ty = imu->halfT * data->g[1];
//    imu->Tz = imu->halfT * data->g[2];

//    imu->kalman->F_t->data[6][7] = imu->kalman->F->data[7][6] = imu->kalman->F_t->data[9][8] = imu->kalman->F->data[8][9] = imu->Tx;
//    imu->kalman->F_t->data[6][8] = imu->kalman->F->data[8][6] = imu->kalman->F_t->data[7][9] = imu->kalman->F->data[9][7] = imu->Ty;
//    imu->kalman->F_t->data[8][7] = imu->kalman->F->data[7][8] = imu->kalman->F_t->data[6][9] = imu->kalman->F->data[9][6] = imu->Tz;

//    imu->kalman->F_t->data[7][6] = imu->kalman->F->data[6][7] = imu->kalman->F_t->data[8][9] = imu->kalman->F->data[9][8] = -imu->Tx;
//    imu->kalman->F_t->data[8][6] = imu->kalman->F->data[6][8] = imu->kalman->F_t->data[9][7] = imu->kalman->F->data[7][9] = -imu->Ty;
//    imu->kalman->F_t->data[9][6] = imu->kalman->F->data[6][9] = imu->kalman->F_t->data[7][8] = imu->kalman->F->data[8][7] = -imu->Tz;


//    // update G
//    imu->kalman->G->data[4][2] = imu->kalman->G->data[2][1] = imu->kalman->G->data[0][0] = (data->dt * data->dt) * 0.5f;
//    imu->kalman->G->data[5][2] = imu->kalman->G->data[3][1] = imu->kalman->G->data[1][0] = data->dt;

//    kalmanPredict(imu->kalman);

//    /*
//     * ***************************************************
//     * Update Step 2
//     * ***************************************************
//     */

//    // first type of update --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

////    //    // normalize accel
////        imu->recipNorm = invSqrt(data->a[0] * data->a[0] + data->a[1] * data->a[1] + data->a[2] * data->a[2]);
////        data->a[0] *= imu->recipNorm;
////        data->a[1] *= imu->recipNorm;
////        data->a[2] *= imu->recipNorm;

////        // normalize gravitation
////        imu->recipNorm = invSqrt(imu->grav[0] * imu->grav[0] + imu->grav[1] * imu->grav[1] + imu->grav[2] * imu->grav[2]);
////        imu->grav_norm[0] = imu->grav[0] * imu->recipNorm;
////        imu->grav_norm[1] = imu->grav[1] * imu->recipNorm;
////        imu->grav_norm[2] = imu->grav[2] * imu->recipNorm;

////        vecMulCrossAngle(data->a, imu->grav_norm, imu->rotateAxisErr, &imu->angleErr);             // computation vector n_a = a_measured × a_calculated; and ∆θa = acos(a_measured * a_calculated)
////        Quaternion_fromAxisAngle(imu->rotateAxisErr, (imu->angleErr * imu->const_u), &imu->q_ae);  // computation error quaternion q_ae
////        Quaternion_multiply_to_arrayLN(&imu->q_ae, &imu->q_a, imu->kalman->Z->data);
////        kalmanUpdate(imu->kalman);

////    //    // second type of update --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    float e0[3] = {0, };
//    e0[0] = accelCalcPitch(data->a[0],data->a[1],data->a[2]);
//    e0[1] = accelCalcRoll(data->a[0],data->a[1],data->a[2]);
//    e0[2] = 0;
//    Quaternion_fromEulerZYX(e0, &imu->q_ae);
//    imu->kalman->Z->data[0][0] = imu->q_ae.w;
//    imu->kalman->Z->data[1][0] = imu->q_ae.v[0];
//    imu->kalman->Z->data[2][0] = imu->q_ae.v[1];
//    imu->kalman->Z->data[3][0] = imu->q_ae.v[2];
//    kalmanUpdate(imu->kalman);
//    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//    computeGravVector(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], imu->grav, imu->calibrationGrav); // a_calculated
//    Quaternion_set(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], &imu->q_a); // qk
//    //Quaternion_set(imu->kalman->Z->data[0][0], imu->kalman->Z->data[1][0], imu->kalman->Z->data[2][0], imu->kalman->Z->data[3][0], &imu->RES); // qk
//    return KALMAN_OK;
//}



////////*************************************DIRECT QUAT************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data) // only quaternion generate from gyroscope
//{

//    imu->halfT = 0.5f * data->dt_sec;
//    imu->Tx = imu->halfT * data->g[0];
//    imu->Ty = imu->halfT * data->g[1];
//    imu->Tz = imu->halfT * data->g[2];

////    imu->quat[0] = imu->quat[0] - (imu->Tx * imu->quat[1]) - (imu->Ty * imu->quat[2]) - (imu->Tz * imu->quat[3]);
////    imu->quat[1] = imu->quat[1] + (imu->Tx * imu->quat[0]) + (imu->Tz * imu->quat[2]) - (imu->Ty * imu->quat[3]);
////    imu->quat[2] = imu->quat[2] + (imu->Ty * imu->quat[0]) - (imu->Tz * imu->quat[1]) + (imu->Tx * imu->quat[3]);
////    imu->quat[3] = imu->quat[3] + (imu->Tz * imu->quat[0]) + (imu->Ty * imu->quat[1]) - (imu->Tx * imu->quat[2]);
////    Quaternion_set(imu->quat[0], imu->quat[1], imu->quat[2], imu->quat[3], &imu->q_a);



//    Quaternion_set(0.0f, data->g[0], data->g[1], data->g[2], &imu->q_ae);
//    Quaternion_multiply(&imu->q_a, &imu->q_ae, &imu->q_i);
//    Quaternion_scalar_multiplication(&imu->q_i, imu->halfT, &imu->q_i);
//    Quaternion_add(&imu->q_a, &imu->q_i, &imu->q_a);
//    Quaternion_normalize(&imu->q_a, &imu->q_a);
//    return KALMAN_OK;
//}


//**************************INVERSE QUAT***********************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data) // only quaternion generate from gyroscope
//{

//    //    imu->halfT = -0.5f * data->dt;
//    //    imu->Tx = imu->halfT * data->g[0];
//    //    imu->Ty = imu->halfT * data->g[1];
//    //    imu->Tz = imu->halfT * data->g[2];

//    //    imu->quat[0] = imu->quat[0] - (imu->Tx * imu->quat[1]) - (imu->Ty * imu->quat[2]) - (imu->Tz * imu->quat[3]);
//    //    imu->quat[1] = imu->quat[1] + (imu->Tx * imu->quat[0]) - (imu->Tz * imu->quat[2]) + (imu->Ty * imu->quat[3]);
//    //    imu->quat[2] = imu->quat[2] + (imu->Ty * imu->quat[0]) + (imu->Tz * imu->quat[1]) - (imu->Tx * imu->quat[3]);
//    //    imu->quat[3] = imu->quat[3] + (imu->Tz * imu->quat[0]) - (imu->Ty * imu->quat[1]) + (imu->Tx * imu->quat[2]);
//    //    Quaternion_set(imu->quat[0], -imu->quat[1], -imu->quat[2], -imu->quat[3], &imu->q_a);

//    imu->halfT = 0.5f * data->dt;
//    imu->Tx = imu->halfT * data->g[0];
//    imu->Ty = imu->halfT * data->g[1];
//    imu->Tz = imu->halfT * data->g[2];
//    imu->quat[0] = imu->quat[0] + (imu->Tx * imu->quat[1]) + (imu->Ty * imu->quat[2]) + (imu->Tz * imu->quat[3]);
//    imu->quat[1] = imu->quat[1] - (imu->Tx * imu->quat[0]) + (imu->Tz * imu->quat[2]) - (imu->Ty * imu->quat[3]);
//    imu->quat[2] = imu->quat[2] - (imu->Ty * imu->quat[0]) - (imu->Tz * imu->quat[1]) + (imu->Tx * imu->quat[3]);
//    imu->quat[3] = imu->quat[3] - (imu->Tz * imu->quat[0]) + (imu->Ty * imu->quat[1]) - (imu->Tx * imu->quat[2]);
//    Quaternion_set(imu->quat[0], -imu->quat[1], -imu->quat[2], -imu->quat[3], &imu->q_a);


//    //    Quaternion_conjugate(&imu->q_a, &imu->q_a);

//    //    Quaternion_set(0.0f, data->g[0], data->g[1], data->g[2], &imu->q_ae);
//    //    Quaternion_multiply(&imu->q_ae, &imu->q_a, &imu->RES);
//    //    Quaternion_scalar_multiplication(&imu->RES, imu->halfT, &imu->RES);
//    //    Quaternion_add(&imu->q_a, &imu->RES, &imu->q_a);

//    //    Quaternion_conjugate(&imu->q_a, &imu->q_a);
//    return KALMAN_OK;
//}

////////*************************************************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data)  // only quaternion generate from accelerometer (yaw pith raw)
//{
//    float e0[3] = {0, };
//    e0[0] = accelCalcPitch(data->a[0],data->a[1],data->a[2]);
//    e0[1] = accelCalcRoll(data->a[0],data->a[1],data->a[2]);
//    e0[2] = 0;
//    Quaternion_fromEulerZYX(e0, &imu->q_a);
//    Quaternion_normalize(&imu->q_a, &imu->q_a);
//    return KALMAN_OK;
//}


#undef KALx
#undef KALz
#undef KALu
