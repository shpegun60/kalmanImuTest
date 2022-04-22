TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    kalman_filter/imu.c \
    kalman_filter/kalman_filter.c \
    main.cpp \
    matrix/matrix.c \
    matrix/matrix_test.c \
    quaternion/quaternion.c \
    quaternion/test_quaternion.c \
    smart_assert/smart_assert.c

HEADERS += \
    kalman_filter/imu.h \
    kalman_filter/kalman_filter.h \
    kalman_filter/kalman_port.h \
    matrix/matrix.h \
    matrix/matrix_test.h \
    quaternion/quaternion.h \
    quaternion/test_quaternion.h \
    smart_assert/smart_assert.h

INCLUDEPATH += \
    kalman_filter\
    quaternion\
    matrix/
    smart_assert/
