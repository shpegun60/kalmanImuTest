/**
 * @file    quaternion.c
 * @brief   A basic quaternion library written in C
 * @date    2019-11-28
 */
#include "quaternion.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>

//-------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

#ifdef FAST_CALCULATE_INV_SQRT
static_assert (sizeof(long) == 4, "need rewrite platform depend invSqrt function or commit FAST_CALCULATE_INV_SQRT because type long != 4");
#endif

float invSqrt(const float x) {
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long *)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float *)&i;
    y = y * (1.5f - (halfx * y * y));
    y = y * (1.5f - (halfx * y * y));
    return y;
}

float fastSqrt(const float n)
{
    /*
     * *********************************************************************************************
     * 1) HIGH precision for math.
     *  -Advantages:
     *  * The convergence of the function is high enough to be used in applications other than games.
     *
     *  -Disadvantages:
     *  * Computing times are much larger.
     *  * The square root of “0” is a number very close to “0” but never “0”.
     *  * The division operation is the bottleneck in this function. because the division operation is more expensive than the other arithmetic operations of Arithmetic Logic Units (ALU)
     **********************************************************************************************
    */
    union {int i; float f;} u;
    u.i = 0x1FB5AD00 + (*(int*)&n >> 1);
    u.f = n / u.f + u.f;
    return n / u.f + u.f * 0.25f;

    /*
     ********************************************************************************************
     * 2) LOW precision for games. (less then origin square)
     * -Advantages:
     *      * The convergence of the function is greater than that of the first method.
     *      * Generates times equal to or greater than the first method.
     *
     * -Disadvantages:
     *      * The square root of “0” is a number very close to “0” but never “0”.
     *******************************************************************************************
    */
    //    union{int i; float f;} u;
    //    u.i = 0x5F375A86 - (*(int*)&n >> 1);
    //    return ((int)3 - n * u.f * u.f) * n * u.f * 0.5f;

    /*
     * ******************************************************************************************
     * 3) LOW precision for games. (more than origin square)
     *  -Advantages:
     *      * When Root of 0 is calculated the function returns 0.
     *      * The convergence of the function is acceptable enough for games.
     *      * It generates very good times.
     *      * The Reciprocal of the root can be calculated by removing the second “n” from the third line. According to the property of: 1 / sqrt (n) * n = sqrt (n).
     *
     *  -Disadvantages:
     *      * Convergence decreases when the root to be calculated is very large.
     * *****************************************************************************************
     */
    //    union {int i; float f;} u;
    //    u.i = 0x2035AD0C + (*(int*)&n >> 1);
    //    return n / u.f + u.f * 0.25f;
}

//-------------------------------------------------------------------------------------------


void Quaternion_set(float w, float v1, float v2, float v3, Quaternion* output)
{
    assert(output != NULL);
    output->w = w;
    output->v[0] = v1;
    output->v[1] = v2;
    output->v[2] = v3;
}

void Quaternion_setIdentity(Quaternion* q)
{
    assert(q != NULL);
    Quaternion_set(1, 0, 0, 0, q);
}

void Quaternion_copy(Quaternion* q, Quaternion* output)
{
    Quaternion_set(q->w, q->v[0], q->v[1], q->v[2], output);
}

bool Quaternion_equal(Quaternion* q1, Quaternion* q2)
{
    bool equalW  = fabs(q1->w - q2->w) <= QUATERNION_EPS;
    bool equalV0 = fabs(q1->v[0] - q2->v[0]) <= QUATERNION_EPS;
    bool equalV1 = fabs(q1->v[1] - q2->v[1]) <= QUATERNION_EPS;
    bool equalV2 = fabs(q1->v[2] - q2->v[2]) <= QUATERNION_EPS;
    return equalW && equalV0 && equalV1 && equalV2;
}

void Quaternion_fprint(FILE* file, Quaternion* q)
{
    fprintf(file, "(%.3f, %.3f, %.3f, %.3f)",
            q->w, q->v[0], q->v[1], q->v[2]);
}


void Quaternion_fromAxisAngle(float axis[3], float angle, Quaternion* output)
{
    assert(output != NULL);
    // Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/
    output->w = cos(angle / 2.0);
    float c = sin(angle / 2.0);
    output->v[0] = c * axis[0];
    output->v[1] = c * axis[1];
    output->v[2] = c * axis[2];
}


float Quaternion_toAxisAngle(Quaternion* q, float output[3])
{
    assert(output != NULL);
    // Formula from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
    float angle = 2.0 * acos(q->w);

    if(q->w < 1.0) {
        // Calculate the axis
#ifdef FAST_CALCULATE_INV_SQRT

        float multiplier = invSqrt(1.0 - q->w * q->w);
#else
        float multiplier = 1.0/sqrt(1.0 - q->w * q->w);
#endif
        output[0] = q->v[0] * multiplier;
        output[1] = q->v[1] * multiplier;
        output[2] = q->v[2] * multiplier;
    } else {
        // Arbitrary normalized axis
        output[0] = 1;
        output[1] = 0;
        output[2] = 0;
    }
    return angle;
}

void Quaternion_fromXRotation(float angle, Quaternion* output)
{
    assert(output != NULL);
    float axis[3] = {1.0, 0, 0};
    Quaternion_fromAxisAngle(axis, angle, output);
}

void Quaternion_fromYRotation(float angle, Quaternion* output)
{
    assert(output != NULL);
    float axis[3] = {0, 1.0, 0};
    Quaternion_fromAxisAngle(axis, angle, output);
}

void Quaternion_fromZRotation(float angle, Quaternion* output)
{
    assert(output != NULL);
    float axis[3] = {0, 0, 1.0};
    Quaternion_fromAxisAngle(axis, angle, output);
}

void Quaternion_fromEulerZYX(float eulerZYX[3], Quaternion* output)
{
    assert(output != NULL);
    // Based on https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    float cy = cos(eulerZYX[2] * 0.5);
    float sy = sin(eulerZYX[2] * 0.5);
    float cr = cos(eulerZYX[0] * 0.5);
    float sr = sin(eulerZYX[0] * 0.5);
    float cp = cos(eulerZYX[1] * 0.5);
    float sp = sin(eulerZYX[1] * 0.5);

    output->w = cy * cr * cp + sy * sr * sp;
    output->v[0] = cy * sr * cp - sy * cr * sp;
    output->v[1] = cy * cr * sp + sy * sr * cp;
    output->v[2] = sy * cr * cp - cy * sr * sp;
}

void Quaternion_toEulerZYX(Quaternion* q, float output[3])
{
    assert(output != NULL);

    // Roll (x-axis rotation)
    float sinr_cosp = +2.0 * (q->w * q->v[0] + q->v[1] * q->v[2]);
    float cosr_cosp = +1.0 - 2.0 * (q->v[0] * q->v[0] + q->v[1] * q->v[1]);
    output[0] = atan2(sinr_cosp, cosr_cosp);

    // Pitch (y-axis rotation)
    float sinp = +2.0 * (q->w * q->v[1] - q->v[2] * q->v[0]);
    if (fabs(sinp) >= 1) {
        output[1] = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    } else {
        output[1] = asin(sinp);
    }

    // Yaw (z-axis rotation)
    float siny_cosp = +2.0 * (q->w * q->v[2] + q->v[0] * q->v[1]);
    float cosy_cosp = +1.0 - 2.0 * (q->v[1] * q->v[1] + q->v[2] * q->v[2]);
    output[2] = atan2(siny_cosp, cosy_cosp);
}

void Quaternion_conjugate(Quaternion* q, Quaternion* output)
{
    assert(output != NULL);
    output->w = q->w;
    output->v[0] = -q->v[0];
    output->v[1] = -q->v[1];
    output->v[2] = -q->v[2];
}

float Quaternion_norm(Quaternion* q)
{
    assert(q != NULL);
#ifdef FAST_CALCULATE_INV_SQRT
    return fastSqrt(q->w*q->w + q->v[0]*q->v[0] + q->v[1]*q->v[1] + q->v[2]*q->v[2]);
#else
    return sqrt(q->w*q->w + q->v[0]*q->v[0] + q->v[1]*q->v[1] + q->v[2]*q->v[2]);
#endif
}

void Quaternion_normalize(Quaternion* q, Quaternion* output)
{
    assert(output != NULL);
#ifdef FAST_CALCULATE_INV_SQRT
    float inverse_len = invSqrt(q->w*q->w + q->v[0]*q->v[0] + q->v[1]*q->v[1] + q->v[2]*q->v[2]);
#else
    float inverse_len = 1.0 / Quaternion_norm(q);
#endif
    Quaternion_set(
                q->w * inverse_len,
                q->v[0] * inverse_len,
            q->v[1] * inverse_len,
            q->v[2] * inverse_len,
            output);
}

void Quaternion_multiply(Quaternion* q1, Quaternion* q2, Quaternion* output)
{
    assert(output != NULL);
    Quaternion result;

    /*
    Formula from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/index.htm
             a*e - b*f - c*g - d*h
        + i (b*e + a*f + c*h- d*g)
        + j (a*g - b*h + c*e + d*f)
        + k (a*h + b*g - c*f + d*e)
    */
    result.w =    q1->w   *q2->w    - q1->v[0]*q2->v[0] - q1->v[1]*q2->v[1] - q1->v[2]*q2->v[2];
    result.v[0] = q1->v[0]*q2->w    + q1->w   *q2->v[0] + q1->v[1]*q2->v[2] - q1->v[2]*q2->v[1];
    result.v[1] = q1->w   *q2->v[1] - q1->v[0]*q2->v[2] + q1->v[1]*q2->w    + q1->v[2]*q2->v[0];
    result.v[2] = q1->w   *q2->v[2] + q1->v[0]*q2->v[1] - q1->v[1]*q2->v[0] + q1->v[2]*q2->w   ;

    *output = result;
}

void Quaternion_multiply_to_array(Quaternion* q1, Quaternion* q2, float* output)
{
    assert(output != NULL);

    /*
    Formula from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/index.htm
             a*e - b*f - c*g - d*h
        + i (b*e + a*f + c*h- d*g)
        + j (a*g - b*h + c*e + d*f)
        + k (a*h + b*g - c*f + d*e)
    */
    output[0] =    q1->w   *q2->w    - q1->v[0]*q2->v[0] - q1->v[1]*q2->v[1] - q1->v[2]*q2->v[2];
    output[1] = q1->v[0]*q2->w    + q1->w   *q2->v[0] + q1->v[1]*q2->v[2] - q1->v[2]*q2->v[1];
    output[2] = q1->w   *q2->v[1] - q1->v[0]*q2->v[2] + q1->v[1]*q2->w    + q1->v[2]*q2->v[0];
    output[3] = q1->w   *q2->v[2] + q1->v[0]*q2->v[1] - q1->v[1]*q2->v[0] + q1->v[2]*q2->w   ;
}

void Quaternion_rotate(Quaternion* q, float v[3], float output[3])
{
    assert(output != NULL);
    float result[3];

    float ww = q->w * q->w;
    float xx = q->v[0] * q->v[0];
    float yy = q->v[1] * q->v[1];
    float zz = q->v[2] * q->v[2];
    float wx = q->w * q->v[0];
    float wy = q->w * q->v[1];
    float wz = q->w * q->v[2];
    float xy = q->v[0] * q->v[1];
    float xz = q->v[0] * q->v[2];
    float yz = q->v[1] * q->v[2];

    // Formula from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm
    // p2.x = w*w*p1.x + 2*y*w*p1.z - 2*z*w*p1.y + x*x*p1.x + 2*y*x*p1.y + 2*z*x*p1.z - z*z*p1.x - y*y*p1.x;
    // p2.y = 2*x*y*p1.x + y*y*p1.y + 2*z*y*p1.z + 2*w*z*p1.x - z*z*p1.y + w*w*p1.y - 2*x*w*p1.z - x*x*p1.y;
    // p2.z = 2*x*z*p1.x + 2*y*z*p1.y + z*z*p1.z - 2*w*y*p1.x - y*y*p1.z + 2*w*x*p1.y - x*x*p1.z + w*w*p1.z;

    result[0] = ww*v[0] + 2*wy*v[2] - 2*wz*v[1] +
            xx*v[0] + 2*xy*v[1] + 2*xz*v[2] -
            zz*v[0] - yy*v[0];
    result[1] = 2*xy*v[0] + yy*v[1] + 2*yz*v[2] +
            2*wz*v[0] - zz*v[1] + ww*v[1] -
            2*wx*v[2] - xx*v[1];
    result[2] = 2*xz*v[0] + 2*yz*v[1] + zz*v[2] -
            2*wy*v[0] - yy*v[2] + 2*wx*v[1] -
            xx*v[2] + ww*v[2];

    // Copy result to output
    output[0] = result[0];
    output[1] = result[1];
    output[2] = result[2];
}

void Quaternion_slerp(Quaternion* q1, Quaternion* q2, float t, Quaternion* output)
{
    Quaternion result;

    // Based on http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
    float cosHalfTheta = q1->w*q2->w + q1->v[0]*q2->v[0] + q1->v[1]*q2->v[1] + q1->v[2]*q2->v[2];

    // if q1=q2 or qa=-q2 then theta = 0 and we can return qa
    if (fabs(cosHalfTheta) >= 1.0) {
        Quaternion_copy(q1, output);
        return;
    }

    float halfTheta = acos(cosHalfTheta);
#ifdef FAST_CALCULATE_INV_SQRT
    float sinHalfTheta = invSqrt(1.0 - cosHalfTheta*cosHalfTheta);
#else
    float sinHalfTheta = 1.0 / sqrt(1.0 - cosHalfTheta*cosHalfTheta);
#endif
    // If theta = 180 degrees then result is not fully defined
    // We could rotate around any axis normal to q1 or q2
    if (fabs(sinHalfTheta) > QUATERNION_EPS_INV) {
        result.w = (q1->w * 0.5 + q2->w * 0.5);
        result.v[0] = (q1->v[0] * 0.5 + q2->v[0] * 0.5);
        result.v[1] = (q1->v[1] * 0.5 + q2->v[1] * 0.5);
        result.v[2] = (q1->v[2] * 0.5 + q2->v[2] * 0.5);
        *output = result;
        return;
    }

    // Calculate Quaternion
    float ratioA = sin((1 - t) * halfTheta) * sinHalfTheta;
    float ratioB = sin(t * halfTheta) * sinHalfTheta;
    result.w = (q1->w * ratioA + q2->w * ratioB);
    result.v[0] = (q1->v[0] * ratioA + q2->v[0] * ratioB);
    result.v[1] = (q1->v[1] * ratioA + q2->v[1] * ratioB);
    result.v[2] = (q1->v[2] * ratioA + q2->v[2] * ratioB);

    *output = result;
}
