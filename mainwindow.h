#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include "viewControl.h"
#include "paintArea.h"
#include "sketchrender.h"
#include <QScrollArea>
#include <QComboBox>
#include <iostream>
#include <QTextBrowser>
#include <QLabel>

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


class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    friend class PaintLabel;
    friend class ViewControl;
    friend class PaintArea;
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void creatColorComboBox(QComboBox*);

    bool saveFile(QString fileName, PaintArea *area);
    int is_front = true;
    bool is_drawing = false;


protected:
    void changeEvent(QEvent*);
    void closeEvent(QCloseEvent*);
    void moveEvent (QMoveEvent*);

private:
    int inner_mode, outer_mode;
    Ui::MainWindow *ui;
    QTextBrowser *browser;
    PaintLabel *skectchBigPanel, *refineBigPanel, *finePanel, *modelPanel, *selectPanel, *deformPanel, *leftPanel, *rightPanel;
    PaintLabel *wireLabel, *smoothLabel, *textureLabel, *undoLabel, *redoLabel, *saveLabel;
    PaintLabel *coarseSaveLabel, *coarseLoadLabel, *coarseClearLabel;
    PaintArea *coarsePanel;
    QWidget *highlightOuterBorder, *highlightInnerBorder;
    ViewControl *viewControl;
    Viewer *viewer;
    QFrame *frame;
    QLabel *coarseLabel, *fineLabel;
    sketchRender *rendsketch;

    QString curFile;
    vector< vector<QVector3D>> contours;
    vector< vector<int>> contourId;
    vector< vector<QPoint>> contour2d;

    void action_Run();
    void action_LoadFace();
    void action_SaveModel();
    std::string deformpath = std::string(localfolder()+"/models/front/saved.png");
    std::string coarsepath = std::string(localfolder()+"/models/front/frontface.png");

    bool has_rendered_front = false;
};

class PaintLabel : public QLabel{
public:
    friend class ViewControl;
    PaintLabel(QWidget *parent, MainWindow* parentwindow, int m, int w = PAINT_SIZE, int h = PAINT_SIZE){
        labelwidth = w; labelheight = h;
        QImage image = QImage(w, h, QImage::Format_RGB32);
        image.fill(qRgb(255, 255, 255));
        setImage(image);
    }
    void setImage(QImage image){
        setPixmap(QPixmap::fromImage(QImage(image.scaled(labelwidth, labelheight, Qt::KeepAspectRatio)).convertToFormat(QImage::Format_RGB32)));
    }
    MainWindow* parentwindow;

    int mode;
    int is_front = true; // for CONTO_MODE
    int labelwidth, labelheight;

protected:

};


#endif
