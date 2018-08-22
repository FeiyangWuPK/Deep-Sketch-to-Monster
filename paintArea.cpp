#include "paintArea.h"
#include "mainWindow.h"
#include <QPainter>
#include <QPrintDialog>
#include <QPrinter>
#include <iostream>

typedef const pair<pair<QPoint, QPoint>, int> qqipair;

inline bool disPoint(QPoint a, QPoint b, int c){
    return (a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y()) <=(c*3/4.0)*(c*3/4.0);
}

PaintArea::PaintArea(QFrame* frame, MainWindow* pw): QWidget(frame) {
    image = QImage(PAINT_SIZE,PAINT_SIZE,QImage::Format_RGB32);  // Penal size and format
    backColor = qRgb(255,255,255);    //Penal background color
    image.fill(backColor);
    parentWindow = pw;
    penColor = Qt::black;
    penWidth = 3;
}


void PaintArea::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.drawImage(0,0,image);
}

void PaintArea::wheelEvent(QWheelEvent *event){
    penWidth += event->delta() / 120;
    if (penWidth <= 2) penWidth = 2;
    if (penWidth >= 64) penWidth = 64;
    parentWindow->setPenwidthDisplay(penWidth);
}

void PaintArea::mousePressEvent(QMouseEvent *event) {
    if (parentWindow->inner_mode == SKETCH_FINE_MODE) return;
    parentWindow->reinitmodel();
    parentWindow->highlightInnerBorder->resize(PAINT_SIZE+2,PAINT_SIZE+2);
    parentWindow->highlightInnerBorder->move(PAINT_SIZE/2+VIEW_SIZE+29, PAINT_SIZE/2+19);
    parentWindow->is_drawing = true;
    startPoint = event->pos();
    if (event->button() == Qt::LeftButton){
        penColor = Qt::black;
    } else if (event->button() == Qt::RightButton){
        penWidth = 6;
        penColor = Qt::white;
    }
}
void PaintArea::mouseMoveEvent(QMouseEvent *event) {
    if (parentWindow->inner_mode == SKETCH_FINE_MODE) return;
    if (event->buttons()&Qt::LeftButton || event->buttons()&Qt::RightButton){
        endPoint = event->pos();
        paint(image);
    }
}
void PaintArea::mouseReleaseEvent(QMouseEvent *event) {
    if (parentWindow->inner_mode == SKETCH_FINE_MODE) return;
    if (event->button() == Qt::LeftButton || event->button() == Qt::RightButton){
        endPoint = event->pos();
        paint(image);
        parentWindow->run_caffe();
    }
    penWidth = 3;
}


void PaintArea::paint(QImage &theImage) {
    QPainter pp(&theImage);
    QPen pen = QPen();
    pen.setColor(penColor);
    pen.setWidth(penWidth);
    pp.setPen(pen);
    pp.drawLine(startPoint, endPoint);
    startPoint = endPoint;
    update();
}

void PaintArea::setImageSize(int width, int height) {
    QImage newImage(width,height,QImage::Format_RGB32);
    image = newImage;
    update();
}

void PaintArea::setImageColor(QColor color) {
    backColor = color.rgb();
    image.fill(backColor);
    update();
}

bool PaintArea::saveImage(const QString &fileName, const char *fileFormat) {
    return QImage(image).save(fileName, fileFormat);
}


void PaintArea::loadFromDir(){
    using namespace std;
    QString fileName = QFileDialog::getOpenFileName(this,
           tr("Open Coarse Sketch"),QString::fromStdString(localfolder()+"\\models\\saved\\"),
            tr("PNG File(*.png)"));
    if (fileName.isEmpty()) return;
    doClear();
    image = QImage(fileName).scaled(PAINT_SIZE, PAINT_SIZE, Qt::KeepAspectRatio);
    update();
    parentWindow->run_caffe();
}

void PaintArea::saveToDir(){
    using namespace std;
    QString fileName = QFileDialog::getSaveFileName(this,
           tr("Save Coarse Sketch"), QString::fromStdString(localfolder()+"\\models\\saved\\"),
           tr("PNG File(*.png);;All Files (*)"));
    if (fileName.isEmpty()) return;
    saveImage(fileName, "png");
}

QSize PaintArea::getImageSize() {
    return image.size();
}

void PaintArea::doPrint() {
     QPrinter printer(QPrinter::HighResolution);

     QPrintDialog *printDialog = new QPrintDialog(&printer, this);
     if (printDialog->exec() == QDialog::Accepted) {
         QPainter painter(&printer);
         QRect rect = painter.viewport();
         QSize size = image.size();
         size.scale(rect.size(), Qt::KeepAspectRatio);
         painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
         painter.setWindow(image.rect());
         painter.drawImage(0, 0, image);
     }
 }


void PaintArea::doClear() {
    image.fill(backColor);
    update();
}

void PaintArea::doRefresh() {
    paint(image);
    update();
}

void PaintArea::setPenWidth(int width) {
    penWidth = width;
}

void PaintArea::setPenColor(QColor color) {
    penColor = color;
}
