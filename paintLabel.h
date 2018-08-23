#ifndef PAINTLABEL_H
#define PAINTLABEL_H
#include <QLabel>
<<<<<<< HEAD
=======
#include <QWidget>
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb

#define PAINT_SIZE 256
#define VIEW_SIZE 900

class MainWindow;

<<<<<<< HEAD
class PaintLabel : public QLabel{
public:
    friend class ViewControl;
    PaintLabel(QWidget *parent, MainWindow* parentwindow, int m, int w = PAINT_SIZE, int h = PAINT_SIZE);
    void setImage(QImage image);
=======
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
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
    MainWindow* parentwindow;

    int mode;
    int is_front = true; // for CONTO_MODE
    int labelwidth, labelheight;
<<<<<<< HEAD

protected:

=======
>>>>>>> d148d5907f79639c1b348cec98c15561235a33eb
};

#endif // PAINTLABEL_H
