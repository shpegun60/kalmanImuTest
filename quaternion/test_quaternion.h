#ifndef __QUATERNION_TEST
#define __QUATERNION_TEST

void ASSERT_SAME_DOUBLE(char* description, double x, double y);
void ASSERT_TRUE(char* description, bool check);
void ASSERT_FALSE(char* description, bool check);

void testQuaternion_set(void);

void testQuaternion_setIdentity(void);

void testQuaternion_copy(void);

void testQuaternion_conjugate(void);
void testQuaternion_equal(void);

void testQuaternion_norm(void);
void testQuaternion_normalize(void);

void testQuaternion_fromAxisAngle(void);

void testQuaternion_fromXRotation(void);

void testQuaternion_fromYRotation(void);

void testQuaternion_fromZRotation(void);
void testQuaternion_toAxisAngle(void);

void testQuaternion_fromEulerZYX(void);
void testQuaternion_toEulerZYX(void);

void testQuaternion_multiply(void);
void testQuaternion_rotate(void);
void testQuaternion_slerp(void);
int test_quaternion(void);

#endif /* __QUATERNION_TEST */
