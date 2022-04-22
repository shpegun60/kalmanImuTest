/**
 * @file    IMU.h
 * @brief
 * @date
 */

#ifndef __IMU
#define __IMU

#include "kalman_filter.h"

typedef struct
{
    KalmanFilter* kalman;
} IMU10Dof;


#endif /* __IMU */

