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
    mainWindow.cpp \
    deform.cpp \
    sketchrender.cpp \
    paintArea.cpp \
    folder.cpp \
    deform.cpp \
    folder.cpp \
    main.cpp \
    mainwindow.cpp \
    paintArea.cpp \
    sketchrender.cpp

HEADERS += \
    mainWindow.h \
    deform.h \
    sketchrender.h \
    lineqn.h \
    paintArea.h \
    folder.h \
    caffeHeader.h \
    caffeParser.h \
    deform.h \
    folder.h \
    lineqn.h \
    mainwindow.h \
    paintArea.h \
    paintlabel.h \
    sketchrender.h

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

