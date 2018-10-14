#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QApplication>
#include <QDesktopWidget>
#include <QGuiApplication>
#include <QMainWindow>
#include <QLabel>
#include <QtCore>
#include <QScreen>
#include <QScrollArea>
#include <QComboBox>
#include <iostream>
#include <QTextBrowser>
#include <QRect>

#include "paintLabel.h"

#define RETRI_MODE 11
#define SKETC_MODE 12
#define REFIN_SELEC_MODE 13
#define REFIN_DEFOR_MODE 14
#define MODEL_MODE 15

#define WIRE_LABEL 21
#define SMOOTH_LABEL 22
#define TEXTURE_LABEL 23
#define UNDO_LABEL 24
#define REDO_LABEL 25
#define LEFT_LABEL 26
#define RIGHT_LABEL 27
#define SAVE_LABEL 28
#define COARSE_LOAD_LABEL 291
#define COARSE_CLEAR_LABEL 292
#define COARSE_SAVE_LABEL 293

#define SKETCH_MODE 31
#define REFINE_MODE 32
#define VIEW_MODE 33
#define REFINE_SELECT_MODE 34
#define REFINE_DEFORM_MODE 35
#define SKETCH_FINE_MODE 36

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    friend class PaintLabel;
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    int is_front = true;
    bool is_drawing = false;
private:
    Ui::MainWindow *ui;

    int inner_mode, outer_mode;
    QTextBrowser *browser;
    PaintLabel *skectchBigPanel, *refineBigPanel, *finePanel, *modelPanel, *selectPanel, *deformPanel, *leftPanel, *rightPanel;
    PaintLabel *wireLabel, *smoothLabel, *textureLabel, *undoLabel, *redoLabel, *saveLabel;
    PaintLabel *coarseSaveLabel, *coarseLoadLabel, *coarseClearLabel;
    QWidget *highlightOuterBorder, *highlightInnerBorder;
    QFrame *frame;
    QLabel *coarseLabel, *fineLabel;
    QString curFile;

    bool has_rendered_front = false;
};

#endif // MAINWINDOW_H
