#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QApplication>
#include <QDesktopWidget>
#include <QGuiApplication>
#include <QMainWindow>
<<<<<<< HEAD
<<<<<<< HEAD
=======
#include <QLabel>
#include <QtCore>
#include <QScreen>
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
#include <QLabel>
#include <QtCore>
#include <QScreen>
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
#include <QScrollArea>
#include <QComboBox>
#include <iostream>
#include <QTextBrowser>
<<<<<<< HEAD
<<<<<<< HEAD
#include <QLabel>
#include <QDir>
=======
#include <QRect>

#include "paintLabel.h"
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
#include <QRect>

#include "paintLabel.h"
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb

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

<<<<<<< HEAD
<<<<<<< HEAD
#define PAINT_SIZE 256
#define VIEW_SIZE 900
#define VERTEX_SIZE 11510
#define MIN_NORMAL_DEGREE 85
#define MAX_NORMAL_DEGREE 95
#define PI 3.14159265

using namespace std;
=======
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
namespace Ui {
class MainWindow;
}

<<<<<<< HEAD
<<<<<<< HEAD
class PaintLabel;

class MainWindow : public QMainWindow {
=======
class MainWindow : public QMainWindow
{
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
class MainWindow : public QMainWindow
{
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
    Q_OBJECT

public:
    friend class PaintLabel;
<<<<<<< HEAD
<<<<<<< HEAD
    friend class ViewControl;
    friend class PaintArea;
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void creatColorComboBox(QComboBox*);

=======
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
    int is_front = true;
    bool is_drawing = false;
private:
    Ui::MainWindow *ui;

<<<<<<< HEAD
<<<<<<< HEAD
=======
    int inner_mode, outer_mode;
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
    int inner_mode, outer_mode;
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
    QTextBrowser *browser;
    PaintLabel *skectchBigPanel, *refineBigPanel, *finePanel, *modelPanel, *selectPanel, *deformPanel, *leftPanel, *rightPanel;
    PaintLabel *wireLabel, *smoothLabel, *textureLabel, *undoLabel, *redoLabel, *saveLabel;
    PaintLabel *coarseSaveLabel, *coarseLoadLabel, *coarseClearLabel;
<<<<<<< HEAD
<<<<<<< HEAD
    QLabel *coarsePanel;
    QWidget *highlightOuterBorder, *highlightInnerBorder;
    QFrame *frame;
    QLabel *coarseLabel, *fineLabel;
    //sketchRender *rendsketch;

    QString curFile;
    //vector< vector<QVector3D>> contours;
    //vector< vector<int>> contourId;
    //vector< vector<QPoint>> contour2d;

    void action_Run();
    void action_LoadFace();
    void action_SaveModel();
    std::string deformpath = std::string(QDir::currentPath().toStdString()+"/models/front/saved.png");
    std::string coarsepath = std::string(QDir::currentPath().toStdString()+"/models/front/frontface.png");
=======
    QWidget *highlightOuterBorder, *highlightInnerBorder;
    QFrame *frame;
    QLabel *coarseLabel, *fineLabel;
    QString curFile;
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
    QWidget *highlightOuterBorder, *highlightInnerBorder;
    QFrame *frame;
    QLabel *coarseLabel, *fineLabel;
    QString curFile;
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb

    bool has_rendered_front = false;
};

<<<<<<< HEAD
<<<<<<< HEAD

#endif
=======
#endif // MAINWINDOW_H
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
=======
#endif // MAINWINDOW_H
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
