QT += core gui \
        widgets \
        xml \
        opengl \
        printsupport

DEFINES += GLOG_NO_ABBREVIATED_SEVERITIES \
           USE_OPENCV

TARGET = deepsketch2face
TEMPLATE = app
SOURCES += main.cpp \
    mainwindow.cpp \
    paintLabel.cpp

HEADERS += \
    lineqn.h \
    paintlabel.h \
    mainwindow.h

FORMS += \
    mainWindow.ui
CONFIG += console

# caffe
INCLUDEPATH +=  /usr/local/include \
                /opt/X11/include/ \
                /usr/local/include/eigen3/

RC_ICONS = srcs/favicon.ico

DISTFILES += \
    srcs/favicon.ico

