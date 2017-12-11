#-------------------------------------------------
#
# Project created by QtCreator 2017-11-29T07:03:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = calibratiob
TEMPLATE = app


SOURCES += main.cpp\
    regression.cpp \
    calibration.cpp \
    calibwindow.cpp

HEADERS  += \
    regression.h \
    calibration.h \
    calibwindow.h

FORMS    += mainwindow.ui

CONFIG += c++11

win32:
{
    INCLUDEPATH += \
                C:\Eigen \
                #C:\Boost
}
