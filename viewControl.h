#ifndef VIEW_CONTROL_H
#define VIEW_CONTROL_H

#include "mouseGrabber.h"
#include "mainWindow.h"
#include <QWidget>
#include <QMouseEvent>
#include <QPoint>
#include <QImage>
#include <vector>
#include <map>
#include <QFrame>

#define CONV_SIZE 1000
#define MAXDIST_SELECTLINE 200
#define SILID 4

typedef vector<QPoint> LINE;
typedef QVector2D vec2;

class MainWindow;

class ViewControl : public QWidget{
public:
    friend class MainWindow;
    ViewControl(QFrame* frame, Viewer* viewer, MainWindow* window);
    ViewControl():ViewControl(NULL, NULL, NULL){}
    void setImageSize(int width,int height);

    bool isModified() const { return modified; }
    bool saveImage(const QString &fileName, const char *fileFormat);
    bool openImage(const QString &fileName);

    QSize getImageSize();

    int getPenWidth(){return penWidth;}
    void setFrontmode(bool bfront){is_frontmode = bfront;}
    void setPenWidth(int width);
    void setPenColor(QColor color);
    QImage getImage(){return image;}
    QImage getFrontImage(){return frontimg;}
    QImage getSideImage(){return sideimg;}

    void doClear();

    void changeImage(QImage newImage);
    void setMaxBound(const QPoint &);
    void setSkeMaxBound(const QPoint &);
    void setSideSkeMaxBound(const QPoint &);
    void saveCroppedImage();    // for DRAG_MODE
    void saveUncroppedImage();  // for SHIFT_MODE
    void saveSketchImage();

    void addShiftPairs();       // if we previously select a line in drag mode, for shift mode (silhoutte)
    void findShiftPairs();      // find the shift pairs, for shift mode (silhoutte)
    void doShiftPairs();        // do shift pairs, for shift mode (silhoutte)

    vector<int> select_line_index;
    vector<int> shift_map_index; // the correspoding obj raw vertices index, to be mapped to parentViewer->shift_map

//    void selectLine();         // select interest line
    void doShiftLine(bool is_up, int level);
    void doShiftLinearea(bool is_up, int level);

    void doLocalDeform(bool is_up, int level);
    void doLineDeform(QPoint a, QPoint b);

    void selectArea();         // select interest area
    void findShiftArea();
    void doShiftArea(bool is_up, int level);
    vector<long long> shift_candidates; // the vertices on the shift line/area that user draw, with height*i + j, for area use only
    map<int, pair<int, int>> shift_map_index_pos; // corresponding 2d projection position, for area use only


    vector<QPoint> selectedLine;

    void test2Dshape(const std::vector<QVector2D> & shape_vertices_2D);

    void resampleline(LINE& line);
    void setContours(vector< vector<QPoint>>& contour, vector< vector<int>>& contourid);
    void setSideContours(vector< vector<QPoint>>& contour, vector< vector<int>>& contourid);

    int findSideSketchId();
    void snapSideDeformLine(int sid, int sj, int ej, QPoint s, QPoint e);
    void mapSideDeformline();



    void snapLine(LINE& l, QPoint s, QPoint e);
    void autoMouthEye(int sid);
    void mappingDeformline(int sid);
    void mappingSideDeformline(int sid);
    void mappingEraseline(int sid);
    void mappingSideEraseline(int sid);
    void temporalEraseSketch(QPoint& p);
    void temporalEraseSideSketch(QPoint& p);
    void eraseSketch();
    void eraseSideSketch();
    void findshiftmap(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);
    void findsideshiftmap(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);
    void clearSkethes();
    void clearSideSketches();
    void updatecanvas();

    void selectHandle(QPoint& p);
    void eraseHandle(QPoint& p);

    void clearSelection();
    void projectSelection(bool draw); //draw or erase selection area
    void projectAllPoints();
    void projectNeighborSelection();
    vector<pair<int, int>> rawproject;
    vector<int> selectedneighbor;


    bool iseraseline = false;
    bool isdrawline = false;
    bool isdrawwrinkle = false;
    bool isdrawdeformline = false;
    int erasesketchid = -1;
    int sideerasesketchid = -1;

    void paintfrontimg();
    void paintsideimg();
    bool is_frontmode = true;
    bool is_sideedit = false;

protected:
    void paintEvent(QPaintEvent *);
    void mousePressEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);
    void keyPressEvent(QKeyEvent *);
    void paint(QImage&);

private:
    MainWindow* parentWindow;
    Viewer* parentViewer;
    QImage image, frontimg, sideimg;
    QImage sketchimage;
    bool areaSelected = false, lineSelected = false;

    vector<QPoint> drawPoints;
    LINE deformline, eraseline, wrinkleline;
    LINE sidedeformline, sideeraseline, sidewrinkleline;

    int cursketch = -1;
    vector<LINE> wrinkles, sidewrinkles;
    vector<vector<QPoint>> sketches, sidesketches;
    vector<vector<QColor>> sketchescolor, sidesketchescolor;
    vector<vector<int>> sketchesId, sidesketchesId;
    vector<bool> bsketchfixed, bsidesketchfixed;
    vector<vector<vector<QPoint>>> sketchesbuffer;
    vector<pair<int, int>> erasedPoints, sideerasedPoints;
    vector<QColor> erasedColor, sideerasedColor;
    vector<int> segmentid, sidesegmentid;



    int sample_cnt = 0;

    QPoint lineStartPoint, startPoint, endPoint;
    bool modified;

    QColor penColor, skeColor, selColor;
    int penWidth, skeWidth;

    QImage tempImage;
    bool isDrawing;

    //for convolve
    QImage finalImage;
    int leftmax, rightmax, upmax, downmax;
    int skeleftmax, skerightmax, skeupmax, skedownmax;
    int sideskeleftmax, sideskerightmax, sideskeupmax, sideskedownmax;
    bool shifted = false;

    //for shift area
    int area_x_left, area_x_right, area_y_up, area_y_down;

};
#endif
