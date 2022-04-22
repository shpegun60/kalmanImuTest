#include "imu.h"
#include <stdlib.h>
#include "smart_assert/smart_assert.h"


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

IMU10Dof* imuCreate(MAT_TYPE dt, MAT_TYPE* Q10x10, MAT_TYPE* R4x4)
{
    M_Assert_BreakSaveCheck((Q10x10 == NULL || R4x4 == NULL), "imuCreate: not valid input covariance matrics", return NULL);
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
    return imu;
}

int imuProceed(IMU10Dof* imu, IMUinput * data)
{

    /*
     * ***************************************************
     * Predict step 1
     * ***************************************************
     */

    // update F, F_t
    imu->kalman->F_t->data[1][0] = imu->kalman->F->data[0][1] = data->dt;
    imu->kalman->F_t->data[3][2] = imu->kalman->F->data[2][3] = data->dt;
    imu->kalman->F_t->data[5][4] = imu->kalman->F->data[4][5] = data->dt;

    imu->halfT = 0.5f * data->dt;
    imu->Tx = imu->halfT * data->gx;
    imu->Ty = imu->halfT * data->gy;
    imu->Tz = imu->halfT * data->gz;

    imu->kalman->F_t->data[7][6] = imu->kalman->F->data[6][7] = -imu->Tx;
    imu->kalman->F_t->data[8][6] = imu->kalman->F->data[6][8] = -imu->Ty;
    imu->kalman->F_t->data[9][6] = imu->kalman->F->data[6][9] = -imu->Tz;

    imu->kalman->F_t->data[6][7] = imu->kalman->F->data[7][6] = imu->Tx;
    imu->kalman->F_t->data[8][7] = imu->kalman->F->data[7][8] = imu->Tz;
    imu->kalman->F_t->data[9][7] = imu->kalman->F->data[7][9] = -imu->Ty;

    imu->kalman->F_t->data[6][8] = imu->kalman->F->data[8][6] = imu->Ty;
    imu->kalman->F_t->data[7][8] = imu->kalman->F->data[8][7] = -imu->Tz;
    imu->kalman->F_t->data[9][8] = imu->kalman->F->data[8][9] = imu->Tx;

    imu->kalman->F_t->data[6][9] = imu->kalman->F->data[9][6] = imu->Tz;
    imu->kalman->F_t->data[7][9] = imu->kalman->F->data[9][7] = imu->Ty;
    imu->kalman->F_t->data[8][9] = imu->kalman->F->data[9][8] = -imu->Tx;

    // update G
    imu->quadHalfT = (data->dt * data->dt) * 0.5f;
    imu->kalman->G->data[0][0] = imu->quadHalfT;
    imu->kalman->G->data[1][0] = data->dt;
    imu->kalman->G->data[2][1] = imu->quadHalfT;
    imu->kalman->G->data[3][1] = data->dt;
    imu->kalman->G->data[4][2] = imu->quadHalfT;
    imu->kalman->G->data[5][2] = data->dt;

    // update U
    imu->kalman->U->data[0][0] = data->ax - imu->grav[0];
    imu->kalman->U->data[1][0] = data->ay - imu->grav[1];
    imu->kalman->U->data[2][0] = data->az - imu->grav[2];

    kalmanPredict(imu->kalman);

    /*
     * ***************************************************
     * Update Step 2
     * ***************************************************
     */





    return KALMAN_OK;
}






