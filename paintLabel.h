#ifndef PAINTLABEL_H
#define PAINTLABEL_H
#include <QLabel>
#include <QWidget>

#define PAINT_SIZE 256
#define VIEW_SIZE 900

class MainWindow;

class PaintLabel : public QLabel
{
public:
    PaintLabel(QWidget *parent, MainWindow* parentwindow, int m, int w = PAINT_SIZE, int h = PAINT_SIZE): QLabel(parent), parentwindow(parentwindow), mode(m){
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
};

#endif // PAINTLABEL_H
