#include <QGLViewer/qglviewer.h>
#include <QGLViewer/mouseGrabber.h>
#include <QGLViewer/manipulatedFrame.h>
#include <GL/glu.h>
#include <vector>
#include <set>
#include <map>
#include <QtOpenGL>
#include <QGLWidget>
#include <QPoint>
#include "deform.h"
#include "sketchrender.h"


#define PAINT_SIZE 256
#define VIEW_SIZE 900
#define VERTEX_SIZE 11510
#define MIN_NORMAL_DEGREE 85
#define MAX_NORMAL_DEGREE 95
#define PI 3.14159265

typedef QVector2D vec2;
typedef QVector3D vec3;

using namespace std;

class CameraPathPlayer : public qglviewer::MouseGrabber {
public:
  CameraPathPlayer(int nb) : pathNb(nb) {}
  void checkIfGrabsMouse(int x, int y, const qglviewer::Camera* const camera);
  int yPos() { return 25*pathNb; }

protected:
  void mousePressEvent(QMouseEvent* const, qglviewer::Camera* const camera) { camera->playPath(pathNb); }

private:
  int pathNb;
};

class Viewer : public QGLViewer {
public :
    bool loadOBJ(const char * path); // load .obj file
    bool loadfixedId(const char* path);
    bool loadsidefixedId(const char* path);
    bool reloadOBJ(const double* new_vertices);
    bool upsample();
    float distThreshold();
    void upsampleArea(int left, int right, int up, int down, vector<QPoint>& line);

    bool reloadVertices();

    void setMatPrameter(const char* name); // set view mode
    void setEnvironment(unsigned int bgr, unsigned int bgg, unsigned int bgb); // set light environment. three entries correspond to background color;

    // show options
    bool getWire(){return wire;}
    void setWire(bool value){wire=value;}
    bool getSmooth(){return smooth;}
    void setSmooth(bool value){smooth=value;}
    bool getTexture(){return texture;}
    void setTexture(bool value){texture=value;}
    void setMeshInvisible(){ wire = false;  smooth = false ; texture = false; }
    void setMeshVisible(){ wire = false;  smooth = true ; texture = false; }
    QVector3D getCameraPos();

    int currentpos = -1;
    //vertices and normals, "sorted" means sorted by given triagnle vector indices ('f' in .obj)
    vector<vector<int>> vertexIndices; // vertices -> raw_vertices
    vector<int> init_vertexIndices; //inital vertexIndices (w/o upsampling)
    vector<vector<QVector3D>> raw_vertices; //all unique vertices
    vector<vector<QVector3D>> vertices; //sorted vertices (by triangle, may not be unique)
    vector<vector<QVector3D>> raw_point_normals; //all unique vertex normals
    vector<vector<QVector3D>> point_normals; //sorted vertex normals (by triangle, may not be unique)
    vector<vector<QVector3D>> poly_normals; // all sorted polygon (triangle) normals
    vector<vector<set<int>>> neighbors; // all neighbors indices of raw_vertices
    vector<pair<float, float>> vts; //textures
    vector<int> vts_label; //textures

    void resizeVectors();  // limit max currentpos;
    void calcNeighbors();
    void removeAfterCurrentpos(); // remove from currentpos to the end

    // map 2D shape and 3D object
    void calc2Dshape();
    vector<int> shape_vertices_index; //raw index of vertices with MIN_NORMAL_DEGREE ~ MAX_NORMAL_DEGREE normal
    // on screen
    vector<QVector2D> shape_vertices_2D; //2D vertices with MIN_NORMAL_DEGREE ~ MAX_NORMAL_DEGREE normal
    vector<double> shape_vertices_2D_z_value; //deep of 2D vertices with MIN_NORMAL_DEGREE ~ MAX_NORMAL_DEGREE normal
    // in world
    vector<QVector3D> shape_vertices_3D; //3D vertices with MIN_NORMAL_DEGREE ~ MAX_NORMAL_DEGREE normal
    map<int, QVector3D> shift_map; //shape_vertices_index with interest -> new QVector3D. need to use shape_vertices_index to get original raw index
    void deform(bool partial, bool shift, int level, const std::vector<int> & selectedNeighbors); // defrom neighbor/all vertices using shift map, in shift/drag mode
    void handleDeform(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);
    void handleDeform(vector<pair<int, vec2>>& shiftmap, const double* new_vertices, vector<bool>& bfixed);
    void handleDeformAllView(vector<pair<int, vec2>>& frontshiftmap, vector<pair<int, vec2>>& sideshiftmap, const double* new_vertices,
                                     vector<bool>& bfrontfixed, vector<bool>& bsidefixed);

//    void handleDeformSideView(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);
    void handleDeformLeftSide(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);
    void handleDeformRightSide(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed);

    void locateLineArea(int left, int right, int up, int down, vector<QPoint>& line, vector<int>& lpt, vector<float>& lwt);
    void localdeform();
    void localROIDeform();

    void printCurrentView();

    void smoothlinevertices(vector<int>& line);
    void calshiftvts(vector<int>& line, vector<vec3>& vts, float zdir);

    vector<int> fixedid, sidefixid;
    vector<bool> bflag, bdeformflag;

    void smoothSampleArea();
    void smoothSampleAreaNormal();

    //new upsampling and select line
    sketchRender* subd;
    void upsamplelinearea(int left, int right, int up, int down, vector<QPoint>& line);

    //highlighted sorted vertices
    //vector<bool> highlight_vertexIndices;
    vector<bool> hightlight_selection_rawindex;
    vector<bool> hightlight_deformneighbor_rawindex;
    vector<bool> hightlight_init;
    void projectfrontsketches(vector< vector<QVector3D>>& sketches,
                         vector< vector<QPoint>>& lines);
    void projectsidesketches(vector< vector<QVector3D>>& sketches,
                         vector< vector<QPoint>>& lines);
    void projectLeftSideSketches(vector< vector<QVector3D>>& sketches,
                                 vector< vector<QPoint>>& lines);
    void projectRightSideSketches(vector< vector<QVector3D>>& sketches,
                                 vector< vector<QPoint>>& lines);
    void MoveToFront();
    void MoveToSide();
    void MoveToLeftSide();
    void MoveToRightSide();


    void centerize(vector<vec3>& verts);
    void inittransform(vector<vec3>& verts, vector<int>& vid,
                       vector<vec3>& targetvt);

    bool sampled = false;

    bool isROISelection = false;
    bool isHandleSelection = false;

protected :
  virtual void draw();
  virtual void init();
  void displayPlayers();
  void updatePlayers();

private:
  bool wire = false, smooth = true, texture = false;

  CameraPathPlayer** player_;
  int nbPlayers_;
  LapDeform* lapdeform;

  GLfloat fMatAmbientFront[4];
  GLfloat fMatDiffuseFront[4];
  GLfloat fMatSpecularFront[4];
  GLfloat fMatAmbientBack[4];
  GLfloat fMatDiffuseBack[4];
  GLfloat fMatSpecularBack[4];
  GLfloat fMatShininess;
  GLfloat fLightAmbient[4];
  GLfloat fLightDiffuse[4];
  GLfloat fLightSpecular[4];
  GLfloat fLightShininess;
  GLfloat fLightPosition[4];
  GLfloat fMatTransparency;

  QImage textureImg;

  //hard code info
  set<int> backpos;
  const int backnum = 1231;
  const float max_dist = 1.59318;
  const float max_x = 0.743122;
  const float max_y = 1.35017;
  const float max_z = 0.0;
  const float min_x = -0.7605;
  const float min_y = -0.843313;
  const float min_z = -1.21561;
  float current_scale = 1;
};
