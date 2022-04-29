#include "imu.h"
#include <stdlib.h>
#include "smart_assert/smart_assert.h"
#include <assert.h>


#define TO_RAD(x) (((x) / 180.0) * M_PI)
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

IMU10Dof* imuCreate(float const_u, float* gravityConst, MAT_TYPE dt, MAT_TYPE* Q10x10, MAT_TYPE* R4x4)
{
    M_Assert_BreakSaveCheck((gravityConst == NULL), "imuCreate: not valid input gravity const vector", return NULL);
    M_Assert_BreakSaveCheck((Q10x10 == NULL || R4x4 == NULL), "imuCreate: not valid input covariance matrics", return NULL);
    M_Assert_BreakSaveCheck((const_u < 0.0f || const_u > 1.0f), "imuCreate: value of u out of range", return NULL);

    IMU10Dof* imu = (IMU10Dof*)malloc(sizeof(IMU10Dof));
    M_Assert_BreakSaveCheck((imu == NULL), "imuCreate: no memory for allocation structure", return NULL);

    imu->kalman = kalmanCreate(NULL, 10, 4, 3);

    // user init functions ------------------------------------------

    MAT_TYPE F_init [10][10] =
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

    MAT_TYPE X_est_init[10][1] =
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

    MAT_TYPE P_init[10][10] =
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


    MAT_TYPE Tg = ( dt * dt ) / 2.0;
    MAT_TYPE G_init [10][3] =
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


    MAT_TYPE H_init[4][10] =
    {
        {0, 0,  0,  0,  0,  0,  1,  0,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  1,  0,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  1,  0},
        {0, 0,  0,  0,  0,  0,  0,  0,  0,  1}
    };
    matrixInitFromArr(imu->kalman->H, H_init[0]);
    matrixInitFromArr_T(imu->kalman->H_t, H_init[0]);


    matrixInitFromArr(imu->kalman->Q, Q10x10);
    matrixInitFromArr(imu->kalman->R, R4x4);
    imu->const_u = const_u;


    // init else vectors
    for(unsigned i = 0; i < 3; ++i) {
        imu->grav[i] = 0.0f;
        imu->rotateAxisErr[i] = 0.0f;
        imu->grav_norm[i] = 0.0f;

        imu->calibrationGrav[i] = gravityConst[i];
    }

    imu->angleErr = 0.0f;
    Quaternion_setIdentity(&imu->q_ae);
    Quaternion_setIdentity(&imu->q_a);

    imu->Tx = 0.0f;
    imu->Ty = 0.0f;
    imu->Tz = 0.0f;
    imu->halfT = 0.0f;
    imu->recipNorm = 0.0f;

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
    /*
     * ***************************************************
     * Predict step 1
     * ***************************************************
     */

    // update F, F_t
    imu->kalman->F_t->data[5][4] = imu->kalman->F->data[4][5] = imu->kalman->F_t->data[3][2]
            = imu->kalman->F->data[2][3] = imu->kalman->F_t->data[1][0] = imu->kalman->F->data[0][1]
            = data->dt;

    imu->halfT = 0.5f * data->dt;
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
    imu->kalman->G->data[4][2] = imu->kalman->G->data[2][1] = imu->kalman->G->data[0][0] = (data->dt * data->dt) * 0.5f;
    imu->kalman->G->data[5][2] = imu->kalman->G->data[3][1] = imu->kalman->G->data[1][0] = data->dt;

    // update U
    for (unsigned int i  = 0; i < 3; ++i) {
        imu->kalman->U->data[i][0] = data->a[i] - imu->grav[i];
    }

    kalmanPredict(imu->kalman);

    /*
     * ***************************************************
     * Update Step 2
     * ***************************************************
     */

    // first type of update --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // normalize accel
    imu->recipNorm = invSqrt(data->a[0] * data->a[0] + data->a[1] * data->a[1] + data->a[2] * data->a[2]);
    data->a[0] *= imu->recipNorm;
    data->a[1] *= imu->recipNorm;
    data->a[2] *= imu->recipNorm;

    // normalize gravitation
    imu->recipNorm = invSqrt(imu->grav[0] * imu->grav[0] + imu->grav[1] * imu->grav[1] + imu->grav[2] * imu->grav[2]);
    imu->grav_norm[0] = imu->grav[0] * imu->recipNorm;
    imu->grav_norm[1] = imu->grav[1] * imu->recipNorm;
    imu->grav_norm[2] = imu->grav[2] * imu->recipNorm;

    vecMulCrossAngle(data->a, imu->grav_norm, imu->rotateAxisErr, &imu->angleErr);             // computation vector n_a = a_measured × a_calculated; and ∆θa = acos(a_measured * a_calculated)
    Quaternion_fromAxisAngle(imu->rotateAxisErr, (imu->angleErr * imu->const_u), &imu->q_ae);  // computation error quaternion q_ae
    Quaternion_multiply_to_arrayLN(&imu->q_ae, &imu->q_a, imu->kalman->Z->data);
    kalmanUpdate(imu->kalman);

//    // second type of update --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//        float e0[3] = {0, };
//        e0[0] = calcPitch(data->a[0],data->a[1],data->a[2]);
//        e0[1] = calcRoll(data->a[0],data->a[1],data->a[2]);
//        e0[2] = 0;
//        Quaternion_fromEulerZYX(e0, &imu->q_ae);
//        imu->kalman->Z->data[0][0] = imu->q_ae.w;
//        imu->kalman->Z->data[1][0] = imu->q_ae.v[0];
//        imu->kalman->Z->data[2][0] = imu->q_ae.v[1];
//        imu->kalman->Z->data[3][0] = imu->q_ae.v[2];
//        kalmanUpdate(imu->kalman);
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    computeGravVector(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], imu->grav, imu->calibrationGrav); // a_calculated
    Quaternion_set(imu->kalman->X_est->data[6][0], imu->kalman->X_est->data[7][0], imu->kalman->X_est->data[8][0], imu->kalman->X_est->data[9][0], &imu->q_a); // qk
    return KALMAN_OK;
}



////*************************************************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data) // only quaternion generate from gyroscope
//{

//    imu->halfT = 0.5f * data->dt;
//    imu->Tx = imu->halfT * data->g[0];
//    imu->Ty = imu->halfT * data->g[1];
//    imu->Tz = imu->halfT * data->g[2];
//    imu->quat[0] = imu->quat[0] - (imu->Tx * imu->quat[1]) - (imu->Ty * imu->quat[2]) - (imu->Tz * imu->quat[3]);
//    imu->quat[1] = imu->quat[1] + (imu->Tx * imu->quat[0]) + (imu->Tz * imu->quat[2]) - (imu->Ty * imu->quat[3]);
//    imu->quat[2] = imu->quat[2] + (imu->Ty * imu->quat[0]) - (imu->Tz * imu->quat[1]) + (imu->Tx * imu->quat[3]);
//    imu->quat[3] = imu->quat[3] + (imu->Tz * imu->quat[0]) + (imu->Ty * imu->quat[1]) - (imu->Tx * imu->quat[2]);
//    Quaternion_set(imu->quat[0], imu->quat[1], imu->quat[2], imu->quat[3], &imu->q_a);



////    Quaternion_set(0.0f, data->g[0], data->g[1], data->g[2], &imu->q_ae);
////    Quaternion_multiply(&imu->q_a, &imu->q_ae, &imu->RES);
////    Quaternion_scalar_multiplication(&imu->RES, imu->halfT, &imu->RES);
////    Quaternion_add(&imu->q_a, &imu->RES, &imu->q_a);
//    return KALMAN_OK;
//}


//////*************************************************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data)  // only quaternion generate from accelerometer (yaw pith raw)
//{

//    float e0[3] = {0, };
//    e0[0] = calcPitch(data->a[0],data->a[1],data->a[2]);
//    e0[1] = calcRoll(data->a[0],data->a[1],data->a[2]);
//    e0[2] = 0;

//    Quaternion_fromEulerZYX(e0, &imu->q_a);
//    Quaternion_normalize(&imu->q_a, &imu->q_a);
//    return KALMAN_OK;
//}


//////*************************************************************************************************************************************************************************************************************
//int imuProceed(IMU10Dof* imu, IMUinput * data)  // only quaternion generate from accelerometer (yaw pith raw)
//{
//    float e0[3] = {0, };
//    e0[0] = calcPitch(data->a[0],data->a[1],data->a[2]);
//    e0[1] = calcRoll(data->a[0],data->a[1],data->a[2]);
//    e0[2] = 0;
//    Quaternion_fromEulerZYX(e0, &imu->q_a);
//    Quaternion_normalize(&imu->q_a, &imu->q_a);
//    return KALMAN_OK;
//}





