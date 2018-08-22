#ifndef SKETCHRENDER_H
#define SKETCHRENDER_H

#include <QVector3D>
#include <vector>
#include "lineqn.h"

using namespace std;

typedef QVector3D vec4;
typedef QVector3D vec3;
typedef QVector2D vec2;

#define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)
#define PREV(i) ((i)>0 ? (i)-1 : (i)+2)

struct Face {
        int v[3];

        Face() {}
        Face(const int &v0, const int &v1, const int &v2)
            { v[0] = v0; v[1] = v1; v[2] = v2; }
        Face(const int *v_)
            { v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2]; }
        int &operator[] (int i) { return v[i]; }
        const int &operator[] (int i) const { return v[i]; }
        operator const int * () const { return &(v[0]); }
        operator const int * () { return &(v[0]); }
        operator int * () { return &(v[0]); }
        int indexof(int v_) const
        {
            return (v[0] == v_) ? 0 :
                   (v[1] == v_) ? 1 :
                   (v[2] == v_) ? 2 : -1;
        }
    };


class sketchRender{
public:
    sketchRender(){}
    ~sketchRender()
    {
        verts.clear();
        normals.clear();
        faces.clear();
    }
    void initRenderSketch(vector<vec3>& vts,vector<int>& fcs);

    void renderSketch(vec3 viewpos, vector<vec3>& vts, vector<vec3>& nms,
                      vector< vector<vec3>>& sketch, vector< vector<int>>& sketchid, bool bfront);
    void getSketches(vector<vec3>& vts, vector< vector<vec3>>& sketch,
                     vector< vector<int>>& sketchid, int label);

    void upsample(vector<vec3>& vts, vector<int>& fcs, vector<bool>& bflag, vector<bool>& bvertflag, float edge_thr);

private:
    void compute_perview(vec3 viewpos, vector<float> &ndotv, vector<float> &kr,
            vector<float> &sctest_num, vector<float> &sctest_den);

    float find_zero_linear(float val0, float val1);
    float find_zero_hermite(int v0, int v1, float val0, float val1, vec3 &grad0, vec3 &grad1);
    vec3 gradkr(int i, vec3 viewpos);

    void draw_face_isoline2(vec3 viewpos, int v0, int v1, int v2, vector<float> &val,
            vector<float> &test_num, vector<float> &test_den, bool do_hermite, bool do_test);

    void draw_face_isoline(vec3 viewpos, int v0, int v1, int v2, vector<float> &val,
            vector<float> &test_num, vector<float> &test_den, vector<float> &ndotv,
            bool do_bfcull, bool do_hermite, bool do_test);

    void calc_occluding_contours(vec3 viewpos, vector<float> &ndotv, vector<float> &kr);
    void calc_suggestive_contours(vec3 viewpos, vector<float> &ndotv, vector<float> &kr,
            vector<float> &sctest_num, vector<float> &sctest_den);

    void draw_isolines(vec3 viewpos, vector<float> &val, vector<float> &test_num, vector<float> &test_den,
            vector<float> &ndotv, bool do_bfcull, bool do_hermite, bool do_test);

    void rot_coord_sys(const vec3 &old_u, const vec3 &old_v,
                  const vec3 &new_norm, vec3 &new_u, vec3 &new_v);

    void diagonalize_curv(const vec3 &old_u, const vec3 &old_v,
                  float ku, float kuv, float kv, const vec3 &new_norm,
                  vec3 &pdir1, vec3 &pdir2, float &k1, float &k2);
    void proj_curv(const vec3 &old_u, const vec3 &old_v, float old_ku, float old_kuv, float old_kv,
               const vec3 &new_u, const vec3 &new_v, float &new_ku, float &new_kuv, float &new_kv);
    void proj_dcurv(const vec3 &old_u, const vec3 &old_v,
            const vec4 old_dcurv, const vec3 &new_u, const vec3 &new_v, vec4 &new_dcurv);



    void convert_strips();
    void tstrip_build(int f, vector<signed char> &face_avail, vector<int> &todo);
    void need_tstrips();
    void need_pointareas();
    void need_curvatures();
    void need_dcurv();
    void need_acrossedge();
    void need_adjacentfaces();

    void loadKeySketches(const char* frontlabel, const char* sidelabel);
    void loadAllSketches(const char* front, const char* leftside, const char* rightside);


    vec3 loop(int v0, int v1, int v2, int v3);
    vec3 opposite(int f, int v);
    vec3 new_loop_edge(int f1, int f2, int v0, int v1, int v2, int v3);
    vec3 avg_bdy(int v);
    float loop_update_alpha(int n);

    void subdiv(vector<vec3>& vts, vector<Face>& fcs);
    void insert_vert(int f, int e);



private:
    vector<vec3> verts;
    vector<vec3> normals;
    vector<Face> faces;
    float sug_thresh, feature_size;
    vector<vec3> pdir1, pdir2;
    vector<float> curv1, curv2;
    vector<vec4> dcurv;
    vector<vec3> sketch3d;
    vector<int> sketchId;
    vector<int> tstrips;
    vector<float> pointareas;
    vector<vec3> cornerareas;
    vector< vector<int>> adjacentfaces;
    vector< Face> across_edge;
    vector< vector<int>> keyFrontSketches;
    vector< vector<int>> keySideSketches;
    vector< vector<int>> keySketches;
    vector<bool> bfrontValid, bsideValid, bValid;
    vector< vector<int>> leftSideSketches, rightSideSketches;


};

#endif // SKETCHRENDER_H
