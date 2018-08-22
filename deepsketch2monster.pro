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
    mouseGrabber.cpp \
    mainWindow.cpp \
    viewControl.cpp \
    deform.cpp \
    sketchrender.cpp \
    paintArea.cpp \
    folder.cpp

HEADERS += \
    mouseGrabber.h \
    mainWindow.h \
    viewControl.h \
    deform.h \
    sketchrender.h \
    lineqn.h \
    paintArea.h \
    folder.h \
    paintLabel.h

FORMS += \
    mainWindow.ui
CONFIG += console

# caffe
INCLUDEPATH +=  ./tools/caffe-master/src \                      # caffe src
                ./tools/libigl/include \                        # libigl include
                ./tools/libQGLViewer \                          # libQGLViewer include
                /usr/local/include \
                /opt/X11/include/ \
                /usr/local/include/eigen3/

CONFIG(debug, debug|release) {
    LIBS += -lgflagsd \
            -lgflags_nothreadsd \
            -lhdf5_cpp \
            -lhdf5_f90cstub \
            -lhdf5_fortran \
            -lhdf5_hl_cpp \
            -lhdf5_hl_f90cstub \
            -lhdf5_hl_fortran \
            -lhdf5_hl \
            -lhdf5 \
            -lhdf5_tools \
            -lLevelDb \
            -llibglog \
            -llibprotobuf \
            -llmdbD \
            -lopencv_calib3d2410d \
            -lopencv_contrib2410d \
            -lopencv_core2410d \
            -lopencv_features2d2410d \
            -lopencv_flann2410d \
            -lopencv_gpu2410d \
            -lopencv_highgui2410d \
            -lopencv_imgproc2410d \
            -lopencv_legacy2410d \
            -lopencv_ml2410d \
            -lopencv_nonfree2410d \
            -lopencv_objdetect2410d \
            -lopencv_ocl2410d \
            -lopencv_photo2410d \
            -lopencv_stitching2410d \
            -lopencv_superres2410d \
            -lopencv_ts2410d \
            -lopencv_video2410d \
            -lopencv_videostab2410d \
            -lszip \
            -lzlib \
            -lcaffe \
            -lcompute_image_mean \
            -lconvert_imageset \
            -lconvert_mnist_data \

} else {
    LIBS += -lgflags \
            -lgflags_nothreads \
            -lhdf5_hl \
            -lhdf5 \
            -lLevelDb \
            -llibglog \
            -llibprotobuf \
            -llmdb \
            -lopencv_calib3d2410 \
            -lopencv_contrib2410 \
            -lopencv_core2410 \
            -lopencv_features2d2410 \
            -lopencv_flann2410 \
            -lopencv_gpu2410 \
            -lopencv_highgui2410 \
            -lopencv_imgproc2410 \
            -lopencv_legacy2410 \
            -lopencv_ml2410 \
            -lopencv_nonfree2410 \
            -lopencv_objdetect2410 \
            -lopencv_ocl2410 \
            -lopencv_photo2410 \
            -lopencv_stitching2410 \
            -lopencv_superres2410 \
            -lopencv_ts2410 \
            -lopencv_video2410 \
            -lopencv_videostab2410 \
            -lcaffe \
            -lcompute_image_mean \
            -lconvert_imageset \
            -lconvert_mnist_data \

}



RC_ICONS = favicon.ico
LIBS *= -L./tools/libQGLViewer/QGLViewer -framework QGLViewer

