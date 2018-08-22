#ifndef PAINTLABEL_H
#define PAINTLABEL_H

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

class PaintLabel : public QLabel{
public:
    friend class ViewControl;
    PaintLabel(QWidget *parent, MainWindow* parentwindow, int m, int w = PAINT_SIZE, int h = PAINT_SIZE);
    void setImage(QImage image);
    MainWindow* parentwindow;

    int mode;
    int is_front = true; // for CONTO_MODE
    int labelwidth, labelheight;

protected:
    void closeTools();
    void openTools();
    void mouseReleaseEvent(QMouseEvent* event);
};


#endif // PAINTLABEL_H
