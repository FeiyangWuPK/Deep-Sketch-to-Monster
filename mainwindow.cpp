#include <GL/glu.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QtOpenGL>
#include <QGLWidget>
#include <QPainter>
#include <QVariant>
#include <QDesktopWidget>
#include <QTextBrowser>
#include <QFont>
#include <ctime>
#include <QLabel>
#include <QPixmap>
#include <QPalette>
#include <fstream>
#include <iostream>

#include "mainWindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    //caffe_setup();

    // set mainwindow position and size
    setFixedSize(1324,984);
    QDesktopWidget *desktop = QApplication::desktop();
    int screenWidth = desktop->width();
    int screenHeight = desktop->height();
    move((screenWidth/2)-(width()/2),(screenHeight/2)-(height()/2)-45);

    // set label colors
    QPalette blue_palette, black_palette;
    blue_palette.setColor(QPalette::WindowText, Qt::blue);
    black_palette.setColor(QPalette::WindowText, Qt::black);

    // create widgets

    frame = new QFrame;
    setCentralWidget(frame);    //set frontFaceScroll to center

    highlightOuterBorder = new QWidget(frame);
    highlightOuterBorder->resize(652,PAINT_SIZE/4+2);
    highlightOuterBorder->setStyleSheet("background-color: red");
    highlightInnerBorder = new QWidget(frame);
    highlightInnerBorder->resize(PAINT_SIZE/2+2,PAINT_SIZE/2+2);
    highlightInnerBorder->setStyleSheet("background-color: red");

    // up panels
    skectchBigPanel = new PaintLabel(frame, this, SKETCH_MODE, 650, PAINT_SIZE/4);
    skectchBigPanel->setMinimumSize(650, PAINT_SIZE/4);

    skectchBigPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\sketch.png")));
    refineBigPanel = new PaintLabel(frame, this, REFINE_MODE, 650, PAINT_SIZE/4);
    refineBigPanel->setMinimumSize(650, PAINT_SIZE/4);
    refineBigPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\refine.png")));


    // left panels
    modelPanel = new PaintLabel(frame, this, VIEW_MODE, PAINT_SIZE/2, PAINT_SIZE/2);
    modelPanel->setMinimumSize(PAINT_SIZE/2, PAINT_SIZE/2);
    modelPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\viewmodel.png")));

    //right panels
    selectPanel = new PaintLabel(frame, this, REFINE_SELECT_MODE, PAINT_SIZE, PAINT_SIZE/2);
    selectPanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE/2);
    selectPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\select.png")));
    selectPanel->setVisible(false);
    deformPanel = new PaintLabel(frame, this, REFINE_DEFORM_MODE, PAINT_SIZE, PAINT_SIZE/2);
    deformPanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE/2);
    deformPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\deform.png")));
    deformPanel->setVisible(false);
    leftPanel  = new PaintLabel(frame, this, LEFT_LABEL, PAINT_SIZE, PAINT_SIZE/2);
    leftPanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE/2);
    leftPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\left.png")));
    leftPanel->setVisible(false);
    rightPanel  = new PaintLabel(frame, this, RIGHT_LABEL, PAINT_SIZE, PAINT_SIZE/2);
    rightPanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE/2);
    rightPanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\right.png")));
    rightPanel->setVisible(false);


    coarsePanel = new PaintArea(frame, this);
    coarsePanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE);
    finePanel = new PaintLabel(frame, this, SKETCH_FINE_MODE);
    finePanel->setMinimumSize(PAINT_SIZE, PAINT_SIZE);
    finePanel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\fine.png")));
    finePanel->is_front = true;
    coarseLabel = new QLabel(frame);
    coarseLabel->setFont(QFont("Arial", 14, QFont::Normal));
    coarseLabel->setPalette(blue_palette);
    coarseLabel->setText("Coarse");
    fineLabel = new QLabel(frame);
    fineLabel->setFont(QFont("Arial", 14, QFont::Normal));
    fineLabel->setPalette(blue_palette);
    fineLabel->setText("Fine");


    // tool labels
    wireLabel = new PaintLabel(frame, this, WIRE_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    wireLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    wireLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\wire.png")));
    smoothLabel = new PaintLabel(frame, this, SMOOTH_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    smoothLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    smoothLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\smooth.png")));
    textureLabel = new PaintLabel(frame, this, TEXTURE_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    textureLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    textureLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\texture.png")));
    undoLabel = new PaintLabel(frame, this, UNDO_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    undoLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    undoLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\undo.png")));
    redoLabel = new PaintLabel(frame, this, REDO_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    redoLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    redoLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\redo.png")));
    saveLabel  = new PaintLabel(frame, this, SAVE_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    saveLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    saveLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\save.png")));


    coarseLoadLabel = new PaintLabel(frame, this, COARSE_LOAD_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    coarseLoadLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    coarseLoadLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\sketch_load.png")));
    coarseClearLabel = new PaintLabel(frame, this, COARSE_CLEAR_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    coarseClearLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    coarseClearLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\sketch_clear.png")));
    coarseSaveLabel = new PaintLabel(frame, this, COARSE_SAVE_LABEL, PAINT_SIZE/4, PAINT_SIZE/4);
    coarseSaveLabel->setMinimumSize(PAINT_SIZE/4, PAINT_SIZE/4);
    coarseSaveLabel->setImage(QImage(QString::fromStdString(localfolder()+"\\images\\sketch_save.png")));


    // set widget positions


    highlightOuterBorder->move(9, 4);
    highlightInnerBorder->move(9, PAINT_SIZE/4+9);

    skectchBigPanel->move(10, 5);
    refineBigPanel->move(665, 5);


    modelPanel->move(10, PAINT_SIZE/4+10);

    coarseLoadLabel->move(PAINT_SIZE/2+VIEW_SIZE+42, PAINT_SIZE/4+10);
    coarseClearLabel->move(PAINT_SIZE/2+VIEW_SIZE+52+PAINT_SIZE/4+10, PAINT_SIZE/4+10);
    coarseSaveLabel->move(PAINT_SIZE/2+VIEW_SIZE+62+PAINT_SIZE/2+20, PAINT_SIZE/4+10);

    coarsePanel->move(PAINT_SIZE/2+VIEW_SIZE+30, PAINT_SIZE/2+20);
    finePanel->move(PAINT_SIZE/2+VIEW_SIZE+30, PAINT_SIZE/2+PAINT_SIZE+30);
    coarseLabel->move(PAINT_SIZE/2+VIEW_SIZE+35, PAINT_SIZE/2+25);
    fineLabel->move(PAINT_SIZE/2+VIEW_SIZE+35, PAINT_SIZE/2+PAINT_SIZE+35);



    selectPanel->move(PAINT_SIZE/2+VIEW_SIZE+30, PAINT_SIZE/4+10);
    deformPanel->move(PAINT_SIZE/2+VIEW_SIZE+30, 3*PAINT_SIZE/4+20);
    leftPanel->move(PAINT_SIZE/2+VIEW_SIZE+30, 5*PAINT_SIZE/4+30);
    rightPanel->move(PAINT_SIZE/2+VIEW_SIZE+30, 7*PAINT_SIZE/4+40);


    wireLabel->move(PAINT_SIZE/2+254,15+PAINT_SIZE/4);
    textureLabel->move(PAINT_SIZE/2+254+2*PAINT_SIZE/4+20, 15+PAINT_SIZE/4);
    smoothLabel->move(PAINT_SIZE/2+254+PAINT_SIZE/4+10, 15+PAINT_SIZE/4);
    undoLabel->move(PAINT_SIZE/2+254+3*PAINT_SIZE/4+30, 15+PAINT_SIZE/4);
    redoLabel->move(PAINT_SIZE/2+254+4*PAINT_SIZE/4+40, 15+PAINT_SIZE/4);
    saveLabel->move(PAINT_SIZE/2+254+5*PAINT_SIZE/4+50, 15+PAINT_SIZE/4);


    curFile = tr("untitled.png");

    outer_mode = SKETCH_MODE;
    inner_mode = VIEW_MODE;

    rendsketch = new sketchRender;

    cout<<"mainwindow fishied"<<endl;
}


MainWindow::~MainWindow() {
    delete ui;
}


void MainWindow::setPenwidthDisplay(int newsize){
//    if (newsize == -1)
//        ui->action_Size->setIconText(QString::number(currentFace->getPenWidth()));
//    else
//        ui->action_Size->setIconText(QString::number(newsize));;
}

void MainWindow::changeEvent(QEvent *e) {
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

void MainWindow::closeEvent(QCloseEvent*){
    qApp->quit();
}

