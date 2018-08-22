#ifndef PAINT_AREA_H
#define PAINT_AREA_H

#include <QWidget>
#include <QMouseEvent>
#include <QPoint>
#include <QImage>
#include <vector>
#include <QFrame>
#include "folder.h"

#define PAINT true
#define ERASE false

using namespace std;

class MainWindow;

class PaintArea : public QWidget {
public:
    friend class MainWindow;
    PaintArea(QFrame* frame, MainWindow* parentWindow);
    PaintArea():PaintArea(NULL, NULL){}
    void setImageSize(int width,int height);
    void setImageColor(QColor color);

    bool saveImage(const QString &fileName, const char *fileFormat);

    QSize getImageSize();
    void doPrint();

    void doClear();
    void doRefresh();

    int getPenWidth(){return penWidth;}
    void setPenWidth(int width);
    void setPenColor(QColor color);
    void loadFromDir();
    void saveToDir();

    QImage image;



protected:
    void paintEvent(QPaintEvent *);
    void mousePressEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);
    void wheelEvent(QWheelEvent *);

    void paint(QImage& theImage);

private:
    MainWindow* parentWindow;
    QRgb backColor;

    QPoint startPoint,endPoint;

    QColor penColor;
    int penWidth;

};
#endif // PAINT_AREA_H
