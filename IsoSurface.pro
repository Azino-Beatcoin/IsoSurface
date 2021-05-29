#-------------------------------------------------
#
# Project created by QtCreator 2017-03-21T15:15:21
#
#-------------------------------------------------
QT       += core gui opengl

CONFIG+=c++17

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = IsoSurface
TEMPLATE = app

SOURCES += main.cpp \
        mainwindow.cpp drawarea.cpp \
        Triangulation.cpp R3Graph.cpp \
        triFunc.cpp skala.cpp ddd.cpp r2geom.cpp roi.cpp voxelset.cpp \
        voxfunc.cpp mainwindow2.cpp \
        QtLL1Calc/calc.cpp QtLL1Calc/lparser.cpp \
        BinHeap.cpp gauss.cpp \
        exp3ddlg.cpp

HEADERS  += mainwindow.h drawarea.h \
    Triangulation.h R3Graph.h \
    triFunc.h skala.h r2geom.h matrix.h voxelset.h \
    roi.h voxfunc.h \
    QtLL1Calc/calc.h QtLL1Calc/lparser.h \
    BinHeap.h gauss.h \
    exp3ddlg.h

FORMS    += mainwindow.ui exp3ddlg.ui

win32: LIBS += -lOpengl32
