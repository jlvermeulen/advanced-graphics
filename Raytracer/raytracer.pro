#-------------------------------------------------
#
# Project created by QtCreator 2014-11-24T13:57:20
#
#-------------------------------------------------

QT       += core gui\
            opengl

CONFIG   += c++11

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = raytracer
TEMPLATE = app

SOURCES  += main.cpp\
            window.cpp \
            glWidget.cpp \
            ColorD.cpp \
            Intersections.cpp \
            Matrix3x3D.cpp \
            Octree.cpp \
            Vector3D.cpp

HEADERS  += window.h \
            glWidget.h \
            BoundingBox.h \
            ColorD.h \
            Intersections.h \
            Matrix3x3D.h \
            Octree.h \
            RayD.h \
            TriangleD.h \
            Vector3D.h \
            VertexD.h

FORMS    += window.ui
