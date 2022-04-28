INCLUDEPATH += $$PWD
DEPENDPATH += $$PWD	

SOURCES += \
    $$PWD/imu.c \
    $$PWD/kalmanImuTop.cpp \
    $$PWD/kalmanImuTest.cpp \
    $$PWD/kalman_filter.c

HEADERS += \
    $$PWD/imu.h \
    $$PWD/kalmanImuTop.h \
    $$PWD/kalmanImuTest.h \
    $$PWD/kalman_filter.h \
    $$PWD/kalman_port.h
