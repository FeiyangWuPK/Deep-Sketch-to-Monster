#include "paintLabel.h"

PaintLabel::PaintLabel(QWidget *parent, MainWindow* parentwindow, int m, int w, int h){
    labelwidth = w; labelheight = h;
    QImage image = QImage(w, h, QImage::Format_RGB32);
    image.fill(qRgb(255, 255, 255));
    setImage(image);
};

void PaintLabel::setImage(QImage image){
    setPixmap(QPixmap::fromImage(QImage(image.scaled(labelwidth,
                        labelheight, Qt::KeepAspectRatio)).convertToFormat(QImage::Format_RGB32)));
};

