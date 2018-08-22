#include"sketchrender.h"
#include <qmath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "folder.h"

using namespace std;

void sketchRender::loadKeySketches(const char* frontlabel, const char* sidelabel)
{
   bfrontValid.clear();
   bfrontValid.resize(verts.size(), false);
   int nvalid, nsketch;
   ifstream flb(frontlabel);
  /* flb>>nvalid;
   for(int i=0; i<nvalid; i++)
   {
        int id; flb>>id;
        bfrontValid[id] = true;
   }*/

   keyFrontSketches.clear();
   flb>>nsketch;
   for(int i=0; i<nsketch; i++)
   {
        vector<int> sketch;
        int np; flb>>np;
        for(int j=0; j<np; j++)
        {
            int id; flb>>id;
            sketch.push_back(id);
        }
        if(i==8||i==2||i==3)
            std::reverse(sketch.begin(), sketch.end());
        keyFrontSketches.push_back(sketch);
   }
   flb.close();


   bsideValid.clear();
   bsideValid.resize(verts.size(), false);
   ifstream slb(sidelabel);
   slb>>nvalid;
   for(int i=0; i<nvalid; i++)
   {
        int id; slb>>id;
        bsideValid[id] = true;
   }
   keySideSketches.clear();
   slb>>nsketch;
   for(int i=0; i<nsketch; i++)
   {
        vector<int> sketch;
        int np; slb>>np;
        for(int j=0; j<np; j++)
        {
            int id; slb>>id;
            sketch.push_back(id);
        }
        keySideSketches.push_back(sketch);
   }

   slb.close();
}

void sketchRender::loadAllSketches(const char* front, const char* leftside, const char* rightside)
{
    int nsketch;
    ifstream fske(front);
    keyFrontSketches.clear();
    fske>>nsketch;
    for(int i=0; i<nsketch; i++)
    {
         vector<int> sketch;
         int np; fske>>np;
         for(int j=0; j<np; j++)
         {
             int id; fske>>id;
             sketch.push_back(id);
         }
         if(i==8||i==2||i==3)
             std::reverse(sketch.begin(), sketch.end());
         keyFrontSketches.push_back(sketch);
    }
    fske.close();

    ifstream lske(leftside);
    leftSideSketches.clear();
    lske>>nsketch;
    for(int i=0; i<nsketch; i++)
    {
         vector<int> sketch;
         int np; lske>>np;
         for(int j=0; j<np; j++)
         {
             int id; lske>>id;
             sketch.push_back(id);
         }
         leftSideSketches.push_back(sketch);
    }
    lske.close();

    ifstream rske(rightside);
    rightSideSketches.clear();
    rske>>nsketch;
    for(int i=0; i<nsketch; i++)
    {
         vector<int> sketch;
         int np; rske>>np;
         for(int j=0; j<np; j++)
         {
             int id; rske>>id;
             sketch.push_back(id);
         }
         rightSideSketches.push_back(sketch);
    }
    rske.close();
}

void sketchRender::need_adjacentfaces()
{
    int nv = verts.size(), nf = faces.size();

    vector<int> numadjacentfaces(nv);
    for (int i = 0; i < nf; i++) {
        numadjacentfaces[faces[i][0]]++;
        numadjacentfaces[faces[i][1]]++;
        numadjacentfaces[faces[i][2]]++;
    }

    adjacentfaces.resize(verts.size());
    for (int i = 0; i < nv; i++)
        adjacentfaces[i].reserve(numadjacentfaces[i]);

     for (int i = 0; i < nf; i++) {
         for (int j = 0; j < 3; j++)
         adjacentfaces[faces[i][j]].push_back(i);
    }
}

void sketchRender::need_acrossedge()
{
    int nf = faces.size();
    across_edge.resize(nf, Face(-1,-1,-1));

        for (int i = 0; i < nf; i++) {
            for (int j = 0; j < 3; j++) {
                if (across_edge[i][j] != -1)
                    continue;
                int v1 = faces[i][(j+1)%3];
                int v2 = faces[i][(j+2)%3];
                const vector<int> &a1 = adjacentfaces[v1];
                const vector<int> &a2 = adjacentfaces[v2];
                for (int k1 = 0; k1 < a1.size(); k1++) {
                    int other = a1[k1];
                    if (other == i)
                        continue;
                    vector<int>::const_iterator it =
                        find(a2.begin(), a2.end(), other);
                    if (it == a2.end())
                        continue;
                    int ind = (faces[other].indexof(v1)+1)%3;
                    if (faces[other][(ind+1)%3] != v2)
                        continue;
                    across_edge[i][j] = other;
                    across_edge[other][ind] = i;
                    break;
                }
            }
        }
}

void sketchRender::tstrip_build(int f, vector<signed char> &face_avail,
    vector<int> &todo)
{
    Face &v = faces[f];
    if (face_avail[f] == 0) {
        tstrips.push_back(v[0]);
        tstrips.push_back(v[1]);
        tstrips.push_back(v[2]);
        tstrips.push_back(-1);
        face_avail[f] = -1;
        return;
    }

    int score[3];
    for (int i = 0; i < 3; i++) {
        score[i] = 0;
        int ae = across_edge[f][i];
        if (ae == -1 || face_avail[ae] < 0)
            continue;
        score[i]++;
        int next_edge = faces[ae].indexof(v[(i+1)%3]);
        int nae = across_edge[ae][next_edge];
        if (nae == -1 || face_avail[nae] < 0)
            continue;
        score[i]++;
        if (face_avail[ae] == 2)
            score[i]++;
    }

    int best_score = max(max(score[0], score[1]), score[2]);
    int best = (score[0] == best_score) ? 0 :
           (score[1] == best_score) ? 1 : 2;

    int vlast2 = v[ best     ];
    int vlast1 = v[(best+1)%3];
    int vnext  = v[(best+2)%3];
    int dir = 1;
    tstrips.push_back(vlast2);
    tstrips.push_back(vlast1);

    while (1) {
        tstrips.push_back(vnext);
        face_avail[f] = -1;
        for (int j = 0; j < 3; j++) {
            int ae = across_edge[f][j];
            if (ae == -1)
                continue;
            if (face_avail[ae] > 0)
                face_avail[ae]--;
            if (face_avail[ae] == 1)
                todo.push_back(ae);
        }

        f = across_edge[f][faces[f].indexof(vlast2)];
        if (f == -1 || face_avail[f] < 0)
            break;
        vlast2 = vlast1;
        vlast1 = vnext;
        vnext = faces[f][(faces[f].indexof(vlast2)+3+dir)%3];
        dir = -dir;
    }

    tstrips.push_back(-1);
}

void sketchRender::need_tstrips()
{
   int nf = faces.size();

   vector<int> todo;
   vector<signed char> face_avail(nf);
   for (int i = 0; i < nf; i++) {
            face_avail[i] = (across_edge[i][0] != -1) +
                    (across_edge[i][1] != -1) +
                    (across_edge[i][2] != -1);
     if (face_avail[i] == 1)
                todo.push_back(i);
  }

  tstrips.reserve(faces.size() * 2);

  int nstrips = 0;
  int i = 0;
  while (i < faces.size()) {
     int next;
     if (todo.empty()) {
                next = i++;
       } else {
          next = todo.back();
          todo.pop_back();
      }
            if (face_avail[next] < 0)
                continue;
            tstrip_build(next, face_avail, todo);
            nstrips++;
        }
   convert_strips();
}

void sketchRender::convert_strips()
{
    if (tstrips.back() != -1) {
            return;
        }

    int len = 0;
    for (int i = tstrips.size() - 2; i >= 0; i--) {
         if (tstrips[i] == -1) {
            tstrips[i+1] = len;
            len = 0;
            } else {
                    tstrips[i+1] = tstrips[i];
                    len++;
                }
     }
     tstrips[0] = len;
}

void sketchRender::initRenderSketch(vector<vec3>& vts, vector<int>& fcs)
{
    verts = vts;
    faces.clear();
    for(int i=0; i<fcs.size()-1; i+=3)
    {
       Face f;
       f[0] = fcs[i]-1;
       f[1] = fcs[i+1]-1;
       f[2] = fcs[i+2]-1;
       faces.push_back(f);
    }

    need_adjacentfaces();
    need_acrossedge();
    need_tstrips();

    loadAllSketches(std::string(localfolder()+"\\models\\labels\\keysketchfront_label.txt").c_str(),
                    std::string(localfolder()+"\\models\\labels\\leftside_label.txt").c_str(),
                    std::string(localfolder()+"\\models\\labels\\rightside_label.txt").c_str());

}

void sketchRender::renderSketch(vec3 viewpos, vector<vec3>& vts, vector<vec3>& nms,
                                vector< vector<vec3>>& sketch, vector< vector<int>>& sketchid, bool bfront)
{
    vector<float> ndotv, kr;
    vector<float> sctest_num, sctest_den;
    verts = vts; normals = nms;

    if(bfront)
    {
        bValid = bfrontValid;
        keySketches = keyFrontSketches;

        for(int i=0; i<keySketches.size(); i++)
        {
            vector<vec3> ske3d;
            vector<int> skeid;
            for(int j=0; j<keySketches[i].size()-1; j++)
            {
                int v = keySketches[i][j];
                ske3d.push_back(verts[v]);
                skeid.push_back(v);
                v = keySketches[i][j+1];
                ske3d.push_back(verts[v]);
                skeid.push_back(v);
            }

            sketch.push_back(ske3d);
            sketchid.push_back(skeid);
        }
    }
    else
    {
        bValid = bsideValid;
        keySketches = keySideSketches;
        sketch3d.clear();
        sketchId.clear();

        for(int i=0; i<keySketches.size(); i++)
        {
            vector<vec3> ske3d;
            vector<int> skeid;
            for(int j=0; j<keySketches[i].size()-1; j++)
            {
                int v = keySketches[i][j];
                ske3d.push_back(verts[v]);
                skeid.push_back(v);
                v = keySketches[i][j+1];
                ske3d.push_back(verts[v]);
                skeid.push_back(v);
            }

            sketch.push_back(ske3d);
            sketchid.push_back(skeid);
        }

      /*  need_pointareas();
        need_curvatures();
        need_dcurv();

        compute_perview(viewpos, ndotv, kr, sctest_num, sctest_den);
        calc_occluding_contours(viewpos, ndotv, kr);
        calc_suggestive_contours(viewpos, ndotv, kr, sctest_num, sctest_den);

         // append key sketches


        sketch.push_back(sketch3d);
        sketchid.push_back(sketchId);*/

    }
}

void sketchRender::getSketches(vector<vec3>& vts, vector< vector<vec3>>& sketch,
                 vector< vector<int>>& sketchid, int label)
{
    switch(label)
    {
    case 0:
       keySketches = keyFrontSketches; break;
    case 1:
       keySketches = leftSideSketches; break;
    case 2:
       keySketches = rightSideSketches; break;
    }

    //upsample keysketches
    for(int i=0; i<keySketches.size(); i++)
    {
        vector<vec3> ske3d;
        vector<int> skeid;
        for(int j=1; j<keySketches[i].size(); j++)
        {
            int v0, v1;
            v0 = keySketches[i][j-1];
            v1 = keySketches[i][j];
            vec3 a, b, p;
            a = vts[v0]; b = vts[v1];
            p = (a+b)*0.5f;
            ske3d.push_back(a);
            ske3d.push_back(p);
            skeid.push_back(v0);
            skeid.push_back(v0);
        }
        int e = keySketches[i][keySketches[i].size()-1];
        ske3d.push_back(vts[e]);
        skeid.push_back(e);
        sketch.push_back(ske3d);
        sketchid.push_back(skeid);
    }

  /*  for(int i=0; i<keySketches.size(); i++)
    {

        for(int j=0; j<keySketches[i].size(); j++)
        {
            int v = keySketches[i][j];
            ske3d.push_back(vts[v]);
            skeid.push_back(v);
        }

        sketch.push_back(ske3d);
        sketchid.push_back(skeid);
    }*/
}

void sketchRender::compute_perview(vec3 viewpos, vector<float> &ndotv, vector<float> &kr,
        vector<float> &sctest_num, vector<float> &sctest_den)
{
    int nv = verts.size();

    float scthresh = sug_thresh / (feature_size*feature_size);
    bool need_DwKr = true;

    ndotv.resize(nv);
    kr.resize(nv);
    sctest_num.resize(nv);
    sctest_den.resize(nv);

        // Compute quantities at each vertex
    for (int i = 0; i < nv; i++) {
       if(!bValid[i]) continue;

            // Compute n DOT v
       vec3 viewdir = viewpos - verts[i];
       float rlv = 1.0f / viewdir.length();
       viewdir *= rlv;
       ndotv[i] = QVector3D::dotProduct(viewdir, normals[i]);

       float u = QVector3D::dotProduct(viewdir, pdir1[i]), u2 = u*u;
       float v = QVector3D::dotProduct(viewdir, pdir2[i]), v2 = v*v;


            // Note:  this is actually Kr * sin^2 theta
      kr[i] = curv1[i] * u2 + curv2[i] * v2;

      if (!need_DwKr) continue;

            // Use DwKr * sin(theta) / cos(theta) for cutoff test
       sctest_num[i] = u2 * (     u*dcurv[i][0] +
           3.0f*v*dcurv[i][1]) +
           v2 * (3.0f*u*dcurv[i][2] +
           v*dcurv[i][3]);

       float csc2theta = 1.0f / (u2 + v2);

       sctest_num[i] *= csc2theta;
       float tr = (curv2[i] - curv1[i]) *
                u * v * csc2theta;


       sctest_num[i] -= 2.0f * ndotv[i] * tr*tr;

       sctest_den[i] = ndotv[i];
       sctest_num[i] -= scthresh * sctest_den[i];

    }
}

float sketchRender::find_zero_linear(float val0, float val1)
{
    return val0 / (val0 - val1);
}

float sketchRender::find_zero_hermite(int v0, int v1, float val0, float val1, vec3 &grad0, vec3 &grad1)
{
     if (val0 == val1)
            return 0.5f;

     vec3 e01 = verts[v1] - verts[v0];
     float d0 = QVector3D::dotProduct(e01, grad0);
     float d1 = QVector3D::dotProduct(e01, grad1);

     float a = 2 * (val0 - val1) + d0 + d1;
     float b = 3 * (val1 - val0) - 2 * d0 - d1;
     float c = d0, d = val0;

     float sl = 0.0f, sr = 1.0f, valsl = val0, valsr = val1;

     float disc = 4 * b - 12 * a * c;
     if (disc > 0 && a != 0) {
        disc = sqrt(disc);
        float r1 = (-2 * b + disc) / (6 * a);
        float r2 = (-2 * b - disc) / (6 * a);
        if (r1 >= 0 && r1 <= 1 && r2 >= 0 && r2 <= 1) {
            float vr1 = (((a * r1 + b) * r1 + c) * r1) + d;
            float vr2 = (((a * r2 + b) * r2 + c) * r2) + d;
                // When extrema have different signs inside an
                // interval with endpoints with different signs,
                // the middle root is in between the two extrema
            if (vr1 < 0.0f && vr2 >= 0.0f ||
                vr1 > 0.0f && vr2 <= 0.0f) {
                        // 3 roots
                if (r1 < r2) {
                            sl = r1;
                            valsl = vr1;
                            sr = r2;
                            valsr = vr2;
                } else {
                            sl = r2;
                            valsl = vr2;
                            sr = r1;
                            valsr = vr1;
                }
              }
            }
        }

        // Bisection method (constant number of interactions)
        for (int iter = 0; iter < 10; iter++) {
            float sbi = (sl + sr) / 2.0f;
            float valsbi = (((a * sbi + b) * sbi) + c) * sbi + d;

            // Keep the half which has different signs
            if (valsl < 0.0f && valsbi >= 0.0f ||
                valsl > 0.0f && valsbi <= 0.0f) {
                    sr = sbi;
                    valsr = valsbi;
            } else {
                sl = sbi;
                valsl = valsbi;
            }
        }

        return 0.5f * (sl + sr);
}

vec3 sketchRender::gradkr(int i, vec3 viewpos)
{
    vec3 viewdir = viewpos - verts[i];
    float rlen_viewdir = 1.0f / viewdir.length();
    viewdir *= rlen_viewdir;

    float ndotv = QVector3D::dotProduct(viewdir, normals[i]);
    float sintheta = sqrtf(1.0f - sqrtf(ndotv));
    float csctheta = 1.0f / sintheta;
    float u = QVector3D::dotProduct(viewdir, pdir1[i]) * csctheta;
    float v = QVector3D::dotProduct(viewdir, pdir2[i]) * csctheta;
    float kr = curv1[i] * u*u + curv2[i] * v*v;
    float tr = u*v * (curv2[i] - curv1[i]);
    float kt = curv1[i] * (1.0f - u*u) +
               curv2[i] * (1.0f - v*v);
    vec3 w     = u * pdir1[i] + v * pdir2[i];
    vec3 wperp = u * pdir2[i] - v * pdir1[i];
    const vec4 &C = dcurv[i];

    vec3 g = pdir1[i] * (u*u*C[0] + 2.0f*u*v*C[1] + v*v*C[2]) +
             pdir2[i] * (u*u*C[1] + 2.0f*u*v*C[2] + v*v*C[3]) -
            2.0f * csctheta * tr * (rlen_viewdir * wperp +
            ndotv * (tr * w + kt * wperp));
    g *= (1.0f - sqrtf(ndotv));
    g -= 2.0f * kr * sintheta * ndotv * (kr * w + tr * wperp);
    return g;
}

void sketchRender::draw_face_isoline2(vec3 viewpos, int v0, int v1, int v2, vector<float> &val,
        vector<float> &test_num, vector<float> &test_den, bool do_hermite, bool do_test)
{
    // How far along each edge?
        float w10 = do_hermite ?
            find_zero_hermite(v0, v1, val[v0], val[v1],
            gradkr(v0, viewpos), gradkr(v1, viewpos)) :
            find_zero_linear(val[v0], val[v1]);
        float w01 = 1.0f - w10;
        float w20 = do_hermite ?
            find_zero_hermite(v0, v2, val[v0], val[v2],
            gradkr(v0, viewpos), gradkr(v2, viewpos)) :
            find_zero_linear(val[v0], val[v2]);
        float w02 = 1.0f - w20;

      //  std::cout<<w01<<" "<<w10<<std::endl;
        // Points along edges
        vec3 p1 = w01 * verts[v0] + w10 * verts[v1];
        vec3 p2 = w02 * verts[v0] + w20 * verts[v2];

        float test_num1 = 1.0f, test_num2 = 1.0f;
        float test_den1 = 1.0f, test_den2 = 1.0f;
        float z1 = 0.0f, z2 = 0.0f;
        bool valid1 = true;

        if (do_test) {
            // Interpolate to find value of test at p1, p2
            test_num1 = w01 * test_num[v0] + w10 * test_num[v1];
            test_num2 = w02 * test_num[v0] + w20 * test_num[v2];
            if (!test_den.empty()) {
                test_den1 = w01 * test_den[v0] + w10 * test_den[v1];
                test_den2 = w02 * test_den[v0] + w20 * test_den[v2];
            }
            // First point is valid iff num1/den1 is positive,
            // i.e. the num and den have the same sign
            valid1 = ((test_num1 >= 0.0f) == (test_den1 >= 0.0f));
            // There are two possible zero crossings of the test,
            // corresponding to zeros of the num and den
            if ((test_num1 >= 0.0f) != (test_num2 >= 0.0f))
                z1 = test_num1 / (test_num1 - test_num2);
            if ((test_den1 >= 0.0f) != (test_den2 >= 0.0f))
                z2 = test_den1 / (test_den1 - test_den2);
            // Sort and order the zero crossings
            if (z1 == 0.0f)
                z1 = z2, z2 = 0.0f;
            else if (z2 < z1)
                std::swap(z1, z2);
        }

        // If the beginning of the segment was not valid, and
        // no zero crossings, then whole segment invalid
        if (!valid1 && !z1 && !z2)
        {
            return;
        }

        // Draw the valid piece(s)
        int npts = 0;
        if (valid1) {
            sketch3d.push_back(p1);
            if(w01>0.5) sketchId.push_back(v0);
            else sketchId.push_back(v1);
            npts++;
        }
        if (z1) {
            sketch3d.push_back((1.0f - z1) * p1 + z1 * p2);
            if(z1<0.5)
            {
                if(w01>0.5) sketchId.push_back(v0);
                else sketchId.push_back(v1);
            }
            else
            {
                if(w02>0.5) sketchId.push_back(v0);
                else sketchId.push_back(v2);
            }
            npts++;
        }
        if (z2) {
            sketch3d.push_back((1.0f - z2) * p1 + z2 * p2);
            if(z2<0.5)
            {
                if(w01>0.5) sketchId.push_back(v0);
                else sketchId.push_back(v1);
            }
            else
            {
                if(w02>0.5) sketchId.push_back(v0);
                else sketchId.push_back(v2);
            }

            npts++;
        }
        if (npts != 2) {
            sketch3d.push_back(p2);
            if(w02>0.5) sketchId.push_back(v0);
            else sketchId.push_back(v2);
        }
}

void sketchRender::draw_face_isoline(vec3 viewpos, int v0, int v1, int v2, vector<float> &val,
        vector<float> &test_num, vector<float> &test_den, vector<float> &ndotv,
        bool do_bfcull, bool do_hermite, bool do_test)
{
    if(!bValid[v0]||!bValid[v1]||!bValid[v2]) return;
    if(ndotv[v0]>0.5||ndotv[v1]>0.5||ndotv[v2]>0.5) return;

    if (do_bfcull && ndotv[v0] <= 0.0f &&
            ndotv[v1] <= 0.0f && ndotv[v2] <= 0.0f)
            return;

        // Quick reject if derivatives are negative
        if (do_test) {
            if (test_den.empty()) {
                if (test_num[v0] <= 0.0f &&
                    test_num[v1] <= 0.0f &&
                    test_num[v2] <= 0.0f)
                    return;
            } else {
                if (test_num[v0] <= 0.0f && test_den[v0] >= 0.0f &&
                    test_num[v1] <= 0.0f && test_den[v1] >= 0.0f &&
                    test_num[v2] <= 0.0f && test_den[v2] >= 0.0f)
                    return;
                if (test_num[v0] >= 0.0f && test_den[v0] <= 0.0f &&
                    test_num[v1] >= 0.0f && test_den[v1] <= 0.0f &&
                    test_num[v2] >= 0.0f && test_den[v2] <= 0.0f)
                    return;
            }
        }


        // Figure out which val has different sign, and draw
        if (val[v0] < 0.0f && val[v1] >= 0.0f && val[v2] >= 0.0f ||
            val[v0] > 0.0f && val[v1] <= 0.0f && val[v2] <= 0.0f)
            draw_face_isoline2(viewpos, v0, v1, v2, val, test_num, test_den,
            do_hermite, do_test);
        else if (val[v1] < 0.0f && val[v2] >= 0.0f && val[v0] >= 0.0f ||
            val[v1] > 0.0f && val[v2] <= 0.0f && val[v0] <= 0.0f)
            draw_face_isoline2(viewpos, v1, v2, v0, val, test_num, test_den,
            do_hermite, do_test);
        else if (val[v2] < 0.0f && val[v0] >= 0.0f && val[v1] >= 0.0f ||
            val[v2] > 0.0f && val[v0] <= 0.0f && val[v1] <= 0.0f)
            draw_face_isoline2(viewpos, v2, v0, v1, val, test_num, test_den,
            do_hermite, do_test);
}

void sketchRender::calc_occluding_contours(vec3 viewpos, vector<float> &ndotv, vector<float> &kr)
{
    draw_isolines(viewpos, ndotv, kr, vector<float>(), ndotv,
            false, false, true);
}

void sketchRender::calc_suggestive_contours(vec3 viewpos, vector<float> &ndotv, vector<float> &kr,
        vector<float> &sctest_num, vector<float> &sctest_den)
{
    draw_isolines(viewpos, kr, sctest_num, sctest_den, ndotv,
            true, true, true);
}

void sketchRender::draw_isolines(vec3 viewpos, vector<float> &val, vector<float> &test_num, vector<float> &test_den,
        vector<float> &ndotv, bool do_bfcull, bool do_hermite, bool do_test)
{
    const int *t = &tstrips[0];
    const int *stripend = t;
    const int *end = t + tstrips.size();

        // Walk through triangle strips
        while (1) {
            if (t >= stripend){
                if (t >= end)
                    return;
                // New strip: each strip is stored as
                // length followed by indices
                stripend = t + 1 + *t;
                // Skip over length plus first two indices of
                // first face
                t += 3;
            }
            // Draw a line if, among the values in this triangle,
            // at least one is positive and one is negative

            const float &v0 = val[*t], &v1 = val[*(t-1)], &v2 = val[*(t-2)];
            if ((v0 > 0.0f || v1 > 0.0f || v2 > 0.0f) &&
                (v0 < 0.0f || v1 < 0.0f || v2 < 0.0f))
                draw_face_isoline(viewpos, *(t-2), *(t-1), *t,
                val, test_num, test_den, ndotv,
                do_bfcull, do_hermite, do_test);
            t++;
        }
}


void sketchRender::need_pointareas()
{
    int nf = faces.size(), nv = verts.size();

    pointareas.clear();
    pointareas.resize(nv, 0);
    cornerareas.clear();
    cornerareas.resize(nf);

    #pragma omp parallel for
        for (int i = 0; i < nf; i++) {
           // cout<<i<<endl;
            if(!bValid[faces[i][0]]&&!bValid[faces[i][1]]&&!bValid[faces[i][2]])
                continue;
            // Edges
            vec3 e[3] = { verts[faces[i][2]] - verts[faces[i][1]],
                     verts[faces[i][0]] - verts[faces[i][2]],
                     verts[faces[i][1]] - verts[faces[i][0]] };
            // Compute corner weights
            float area = 0.5f * QVector3D::crossProduct(e[0], e[1]).length();
            float l2[3] = { e[0].lengthSquared(), e[1].lengthSquared(), e[2].lengthSquared() };
            float ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
                    l2[1] * (l2[2] + l2[0] - l2[1]),
                    l2[2] * (l2[0] + l2[1] - l2[2]) };
            if (ew[0] <= 0.0f) {
                cornerareas[i][1] = -0.25f * l2[2] * area /
                            QVector3D::dotProduct(e[0], e[2]);
                cornerareas[i][2] = -0.25f * l2[1] * area /
                            QVector3D::dotProduct(e[0], e[1]);
                cornerareas[i][0] = area - cornerareas[i][1] -
                            cornerareas[i][2];
            } else if (ew[1] <= 0.0f) {
                cornerareas[i][2] = -0.25f * l2[0] * area /
                            QVector3D::dotProduct(e[1], e[0]);
                cornerareas[i][0] = -0.25f * l2[2] * area /
                            QVector3D::dotProduct(e[1], e[2]);
                cornerareas[i][1] = area - cornerareas[i][2] -
                            cornerareas[i][0];
            } else if (ew[2] <= 0.0f) {
                cornerareas[i][0] = -0.25f * l2[1] * area /
                            QVector3D::dotProduct(e[2], e[1]);
                cornerareas[i][1] = -0.25f * l2[0] * area /
                            QVector3D::dotProduct(e[2], e[0]);
                cornerareas[i][2] = area - cornerareas[i][0] -
                            cornerareas[i][1];
            } else {
                float ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
                for (int j = 0; j < 3; j++)
                    cornerareas[i][j] = ewscale * (ew[(j+1)%3] +
                                       ew[(j+2)%3]);
            }
    #pragma omp atomic
            pointareas[faces[i][0]] += cornerareas[i][0];
    #pragma omp atomic
            pointareas[faces[i][1]] += cornerareas[i][1];
    #pragma omp atomic
            pointareas[faces[i][2]] += cornerareas[i][2];
        }
}

void sketchRender::rot_coord_sys(const vec3 &old_u, const vec3 &old_v,
              const vec3 &new_norm,
              vec3 &new_u, vec3 &new_v)
{
    new_u = old_u;
    new_v = old_v;
    vec3 old_norm = QVector3D::crossProduct(old_u, old_v);
    float ndot = QVector3D::dotProduct(old_norm, new_norm);
    if (unlikely(ndot <= -1.0f)) {
        new_u = -new_u;
        new_v = -new_v;
        return;
    }
    vec3 perp_old = new_norm - ndot * old_norm;
    vec3 dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
    new_u -= dperp * QVector3D::dotProduct(new_u, perp_old);
    new_v -= dperp * QVector3D::dotProduct(new_v, perp_old);
}

void sketchRender::proj_curv(const vec3 &old_u, const vec3 &old_v,
           float old_ku, float old_kuv, float old_kv,
           const vec3 &new_u, const vec3 &new_v,
           float &new_ku, float &new_kuv, float &new_kv)
{
    vec3 r_new_u, r_new_v;
    rot_coord_sys(new_u, new_v, QVector3D::crossProduct(old_u, old_v), r_new_u, r_new_v);

    float u1 = QVector3D::dotProduct(r_new_u, old_u);
    float v1 = QVector3D::dotProduct(r_new_u, old_v);
    float u2 = QVector3D::dotProduct(r_new_v, old_u);
    float v2 = QVector3D::dotProduct(r_new_v, old_v);

    new_ku  = old_ku * u1*u1 + old_kuv * (2.0f  * u1*v1) + old_kv * v1*v1;
    new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
    new_kv  = old_ku * u2*u2 + old_kuv * (2.0f  * u2*v2) + old_kv * v2*v2;
}


// Like the above, but for dcurv
void sketchRender::proj_dcurv(const vec3 &old_u, const vec3 &old_v,
        const vec4 old_dcurv, const vec3 &new_u, const vec3 &new_v, vec4 &new_dcurv)
{
    vec3 r_new_u, r_new_v;
    rot_coord_sys(new_u, new_v, QVector3D::crossProduct(old_u, old_v), r_new_u, r_new_v);

    float u1 = QVector3D::dotProduct(r_new_u, old_u);
    float v1 = QVector3D::dotProduct(r_new_u, old_v);
    float u2 = QVector3D::dotProduct(r_new_v, old_u);
    float v2 = QVector3D::dotProduct(r_new_v, old_v);

    new_dcurv[0] = old_dcurv[0]*u1*u1*u1 +
               old_dcurv[1]*3.0f*u1*u1*v1 +
               old_dcurv[2]*3.0f*u1*v1*v1 +
               old_dcurv[3]*v1*v1*v1;
    new_dcurv[1] = old_dcurv[0]*u1*u1*u2 +
               old_dcurv[1]*(u1*u1*v2 + 2.0f*u2*u1*v1) +
               old_dcurv[2]*(u2*v1*v1 + 2.0f*u1*v1*v2) +
               old_dcurv[3]*v1*v1*v2;
    new_dcurv[2] = old_dcurv[0]*u1*u2*u2 +
               old_dcurv[1]*(u2*u2*v1 + 2.0f*u1*u2*v2) +
               old_dcurv[2]*(u1*v2*v2 + 2.0f*u2*v2*v1) +
               old_dcurv[3]*v1*v2*v2;
    new_dcurv[3] = old_dcurv[0]*u2*u2*u2 +
               old_dcurv[1]*3.0f*u2*u2*v2 +
               old_dcurv[2]*3.0f*u2*v2*v2 +
               old_dcurv[3]*v2*v2*v2;
}

void sketchRender::diagonalize_curv(const vec3 &old_u, const vec3 &old_v,
              float ku, float kuv, float kv, const vec3 &new_norm,
              vec3 &pdir1, vec3 &pdir2, float &k1, float &k2)
{
    vec3 r_old_u, r_old_v;
    rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

    float c = 1, s = 0, tt = 0;
    if (likely(kuv != 0.0f)) {
        // Jacobi rotation to diagonalize
        float h = 0.5f * (kv - ku) / kuv;
        tt = (h < 0.0f) ?
            1.0f / (h - sqrt(1.0f + h*h)) :
            1.0f / (h + sqrt(1.0f + h*h));
        c = 1.0f / sqrt(1.0f + tt*tt);
        s = tt * c;
    }

    k1 = ku - tt * kuv;
    k2 = kv + tt * kuv;

    if (fabs(k1) >= fabs(k2)) {
        pdir1 = c*r_old_u - s*r_old_v;
    } else {
        std::swap(k1, k2);
        pdir1 = s*r_old_u + c*r_old_v;
    }
    pdir2 = QVector3D::crossProduct(new_norm, pdir1);
}


void sketchRender::need_curvatures()
{
        // Resize the arrays we'll be using
     int nv = verts.size(), nf = faces.size();
     curv1.clear(); curv1.resize(nv); curv2.clear(); curv2.resize(nv);
     pdir1.clear(); pdir1.resize(nv); pdir2.clear(); pdir2.resize(nv);
     vector<float> curv12(nv);

        // Set up an initial coordinate system per vertex
     for (int i = 0; i < nf; i++) {
            if(!bValid[faces[i][0]]&&!bValid[faces[i][1]]&&!bValid[faces[i][2]])
                continue;
            pdir1[faces[i][0]] = verts[faces[i][1]] -
                         verts[faces[i][0]];
            pdir1[faces[i][1]] = verts[faces[i][2]] -
                         verts[faces[i][1]];
            pdir1[faces[i][2]] = verts[faces[i][0]] -
                         verts[faces[i][2]];
        }
    #pragma omp parallel for
        for (int i = 0; i < nv; i++) {
            if(!bValid[i]) continue;
            pdir1[i] =QVector3D::crossProduct(pdir1[i], normals[i]);
            pdir1[i].normalize();
            pdir2[i]=QVector3D::crossProduct(normals[i], pdir1[i]);
        }

        // Compute curvature per-face
    #pragma omp parallel for
        for (int i = 0; i < nf; i++) {
            if(!bValid[faces[i][0]]&&!bValid[faces[i][1]]&&!bValid[faces[i][2]])
                continue;
            // Edges
            vec3 e[3] = { verts[faces[i][2]] - verts[faces[i][1]],
                     verts[faces[i][0]] - verts[faces[i][2]],
                     verts[faces[i][1]] - verts[faces[i][0]] };

            // N-T-B coordinate system per face
            vec3 t = e[0];
            t.normalize();
            vec3 n = QVector3D::crossProduct(e[0], e[1]);
            vec3 b = QVector3D::crossProduct(n, t);
            b.normalize();

            // Estimate curvature based on variation of normals
            // along edges
            float m[3] = { 0, 0, 0 };
            float w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
            for (int j = 0; j < 3; j++) {
                float u = QVector3D::dotProduct(e[j], t);
                float v = QVector3D::dotProduct(e[j], b);
                w[0][0] += u*u;
                w[0][1] += u*v;
                //w[1][1] += v*v + u*u;
                //w[1][2] += u*v;
                w[2][2] += v*v;
                vec3 dn = normals[faces[i][PREV(j)]] -
                     normals[faces[i][NEXT(j)]];
                float dnu = QVector3D::dotProduct(dn, t);
                float dnv = QVector3D::dotProduct(dn, b);
                m[0] += dnu*u;
                m[1] += dnu*v + dnv*u;
                m[2] += dnv*v;
            }
            w[1][1] = w[0][0] + w[2][2];
            w[1][2] = w[0][1];

            // Least squares solution
            float diag[3];
            if (!ldltdc<float,3>(w, diag)) {
                //dprintf("ldltdc failed!\n");
                continue;
            }
            ldltsl<float,3>(w, diag, m, m);

            // Push it back out to the vertices
            for (int j = 0; j < 3; j++) {
                int vj = faces[i][j];
                if(!bValid[vj]) continue;
                float c1, c12, c2;
                proj_curv(t, b, m[0], m[1], m[2],
                      pdir1[vj], pdir2[vj], c1, c12, c2);
                float wt = cornerareas[i][j] / pointareas[vj];
    #pragma omp atomic
                curv1[vj]  += wt * c1;
    #pragma omp atomic
                curv12[vj] += wt * c12;
    #pragma omp atomic
                curv2[vj]  += wt * c2;
            }
        }

        // Compute principal directions and curvatures at each vertex
    #pragma omp parallel for
        for (int i = 0; i < nv; i++) {
            if(!bValid[i]) continue;
            diagonalize_curv(pdir1[i], pdir2[i],
                     curv1[i], curv12[i], curv2[i],
                     normals[i], pdir1[i], pdir2[i],
                     curv1[i], curv2[i]);
        }
}

void sketchRender::need_dcurv()
{
        // Resize the arrays we'll be using
    int nv = verts.size(), nf = faces.size();
        dcurv.clear(); dcurv.resize(nv);

        // Compute dcurv per-face
    #pragma omp parallel for
        for (int i = 0; i < nf; i++) {
            if(!bValid[faces[i][0]]&&!bValid[faces[i][1]]&&!bValid[faces[i][2]])
                continue;
            // Edges
            vec3 e[3] = { verts[faces[i][2]] - verts[faces[i][1]],
                     verts[faces[i][0]] - verts[faces[i][2]],
                     verts[faces[i][1]] - verts[faces[i][0]] };

            // N-T-B coordinate system per face
            vec3 t = e[0];
            t.normalize();
            vec3 n = QVector3D::crossProduct(e[0], e[1]);
            vec3 b = QVector3D::crossProduct(n, t);
            b.normalize();

            // Project curvature tensor from each vertex into this
            // face's coordinate system
            vec3 fcurv[3];
            for (int j = 0; j < 3; j++) {
                int vj = faces[i][j];
                if(!bValid[vj]) continue;
                proj_curv(pdir1[vj], pdir2[vj], curv1[vj], 0, curv2[vj],
                      t, b, fcurv[j][0], fcurv[j][1], fcurv[j][2]);

            }

            // Estimate dcurv based on variation of curvature along edges
            float m[4] = { 0, 0, 0, 0 };
            float w[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };
            for (int j = 0; j < 3; j++) {
                int vj = faces[i][j];
                if(!bValid[vj]) continue;
                // Variation of curvature along each edge
                vec3 dfcurv = fcurv[PREV(j)] - fcurv[NEXT(j)];
                float u = QVector3D::dotProduct(e[j], t);
                float v = QVector3D::dotProduct(e[j], b);
                float u2 = u*u, v2 = v*v, uv = u*v;
                w[0][0] += u2;
                w[0][1] += uv;
                w[3][3] += v2;
                m[0] += u*dfcurv[0];
                m[1] += v*dfcurv[0] + 2.0f*u*dfcurv[1];
                m[2] += 2.0f*v*dfcurv[1] + u*dfcurv[2];
                m[3] += v*dfcurv[2];
            }
            w[1][1] = 2.0f * w[0][0] + w[3][3];
            w[1][2] = 2.0f * w[0][1];
            w[2][2] = w[0][0] + 2.0f * w[3][3];
            w[2][3] = w[0][1];

            // Least squares solution
            float d[4];
            if (!ldltdc<float,4>(w, d)) {
                //dprintf("ldltdc failed!\n");
                continue;
            }
            ldltsl<float,4>(w, d, m, m);
            vec4 face_dcurv;
            face_dcurv[0] = m[0];
            face_dcurv[1] = m[1];
            face_dcurv[2] = m[2];
            face_dcurv[3] = m[3];


            // Push it back out to each vertex
            for (int j = 0; j < 3; j++) {
                int vj = faces[i][j];
                if(!bValid[vj]) continue;
                vec4 this_vert_dcurv;
                proj_dcurv(t, b, face_dcurv,
                       pdir1[vj], pdir2[vj], this_vert_dcurv);
                float wt = cornerareas[i][j] / pointareas[vj];
                dcurv[vj] += wt * this_vert_dcurv;
            }
        }
}

vec3 sketchRender::loop(int v0, int v1, int v2, int v3)
{
    return 0.125f * (verts[v0] + verts[v3]) +
           0.375f * (verts[v1] + verts[v2]);
}

vec3 sketchRender::avg_bdy(int v)
{
    vec3 p(0, 0, 0);
    int n = 0;
    const vector<int> &a = adjacentfaces[v];
    for (size_t i = 0; i < a.size(); i++) {
        int f = a[i];
        for (int j = 0; j < 3; j++) {
            if (across_edge[f][j] == -1) {
                p += verts[faces[f][NEXT(j)]];
                p += verts[faces[f][PREV(j)]];
                n += 2;
            }
        }
    }
    return p * (1.0f / n);
}

float sketchRender::loop_update_alpha(int n)
{
    if (n == 3)	  return 0.3438f;
    else if (n == 4)  return 0.4625f;
    else if (n == 5)  return 0.5625f;
    else              return 0.625f;
}

void sketchRender::insert_vert(int f, int e)
{
    int v1 = faces[f][NEXT(e)], v2 = faces[f][PREV(e)];

        int ae = across_edge[f][e];
        if (ae == -1) {
            // Boundary
            vec3 p = 0.5f * (verts[v1] + verts[v2]);
            verts.push_back(p);
            return;
        }

        int v0 = faces[f][e];
        const Face &aef = faces[ae];
        int v3 = aef[NEXT(aef.indexof(v1))];
        vec3 p;
        p = loop(v0, v1, v2, v3);

        verts.push_back(p);
}

void sketchRender::subdiv(vector<vec3>& vts, vector<Face>& fcs)
{
    verts = vts; faces = fcs;
    adjacentfaces.clear(); across_edge.clear();
    need_adjacentfaces(); need_acrossedge();

    // Introduce new vertices
    int nf = faces.size();
    vector<Face> newverts(nf, Face(-1,-1,-1));
    int old_nv = verts.size();
    verts.reserve(4 * old_nv);
    vector<int> newvert_count(old_nv + 3*nf);

    for (int i = 0; i < nf; i++) {
         for (int j = 0; j < 3; j++) {
           if (newverts[i][j] != -1)
                    continue;
                int ae = across_edge[i][j];
                if (ae != -1) {
                    if (across_edge[ae][0] == i)
                        newverts[i][j] = newverts[ae][0];
                    else if (across_edge[ae][1] == i)
                        newverts[i][j] = newverts[ae][1];
                    else if (across_edge[ae][2] == i)
                        newverts[i][j] = newverts[ae][2];
                }
                if (newverts[i][j] != -1)
                    continue;

                insert_vert(i, j);
                newverts[i][j] = verts.size() - 1;
                if (ae != -1) {
                    if (across_edge[ae][0] == i)
                        newverts[ae][0] = newverts[i][j];
                    else if (across_edge[ae][1] == i)
                        newverts[ae][1] = newverts[i][j];
                    else if (across_edge[ae][2] == i)
                        newverts[ae][2] = newverts[i][j];
                }
            }
        }



        vector<vec3> oldvertices = verts;
        for (int i = 0; i < old_nv; i++) {
                vec3 bdyavg, nbdyavg;
                int nbdy = 0, nnbdy = 0;
                int naf = adjacentfaces[i].size();
                if (!naf)
                    continue;
                for (int j = 0; j < naf; j++) {
                    int af = adjacentfaces[i][j];
                    int afi = faces[af].indexof(i);
                    int n1 = NEXT(afi);
                    int n2 = PREV(afi);
                    if (across_edge[af][n1] == -1) {
                        bdyavg += oldvertices[faces[af][n2]];
                        nbdy++;
                    } else {
                        nbdyavg += oldvertices[faces[af][n2]];
                        nnbdy++;
                    }
                    if (across_edge[af][n2] == -1) {
                        bdyavg += oldvertices[faces[af][n1]];
                        nbdy++;
                    } else {
                        nbdyavg += oldvertices[faces[af][n1]];
                        nnbdy++;
                    }
                }

                float alpha;
                vec3 newpt;
                if (nbdy) {
                    newpt = bdyavg / (float) nbdy;
                    alpha = 0.75f;
                } else if (nnbdy) {
                    newpt = nbdyavg / (float) nnbdy;
                    alpha = loop_update_alpha(nnbdy/2);
                } else {
                    continue;
                }
                verts[i] *= alpha;
                verts[i] += (1.0f - alpha) * newpt;
            }

        // Insert new faces
        adjacentfaces.clear(); across_edge.clear();
        faces.reserve(4*nf);
        for (int i = 0; i < nf; i++) {
            Face &v = faces[i];
            Face &n = newverts[i];
            faces.push_back(Face(v[0], n[2], n[1]));
            faces.push_back(Face(v[1], n[0], n[2]));
            faces.push_back(Face(v[2], n[1], n[0]));
            v = n;
        }

        vts = verts; fcs = faces;
}

void sketchRender::upsample(vector<vec3>& vts, vector<int>& fcs, vector<bool>& bflag, vector<bool>& bvertflag, float edge_thr)
{
    verts = vts;
    faces.clear();
    for(int i=0; i<fcs.size()-1; i+=3)
    {
       Face f;
       f[0] = fcs[i]-1;
       f[1] = fcs[i+1]-1;
       f[2] = fcs[i+2]-1;
       faces.push_back(f);
    }

    adjacentfaces.clear(); across_edge.clear();
    need_adjacentfaces(); need_acrossedge();

    int nf = faces.size();
    vector< vector<int>> evid;
    evid.resize(nf);
    for(int f=0; f<nf; f++)
        evid[f].resize(3, -1);

    vector<Face> tfaces = faces;

    for(int f=0; f<nf; f++)
    {
        int a, b, c;
        a = faces[f][0]; b = faces[f][1]; c = faces[f][2];
        if(!bflag[a]||!bflag[b]||!bflag[c]) continue;

        vec3 va, vb, vc;
        va = vts[a]; vb = vts[b]; vc = vts[c];
        float lab, lbc, lac;
        lab =(va-vb).length();
        lbc = (vb-vc).length();
        lac = (va-vc).length();
        if(lab<edge_thr && lbc<edge_thr && lac<edge_thr) continue;

        vector<int> nvt; nvt.resize(3, -1);
        int af = across_edge[f][0];
        int bf = across_edge[f][1];
        int cf = across_edge[f][2];

        if(af>=0 && lbc>edge_thr)
        {
            int ii = faces[af].indexof(b);
            int jj = faces[af].indexof(c);
            int ov = faces[af][3-ii-jj];
            if(bflag[ov])
            {
                int id = evid[f][0];
                if(id>=0)
                    nvt[0] = id;
                else
                {
                  //  verts.push_back((vb+vc)*0.5f);
                    verts.push_back((vc+vb)*0.375f+(verts[a]+verts[ov])*0.125f);
                    bflag.push_back(true); bvertflag.push_back(true);
                    nvt[0] = verts.size()-1;
                    evid[f][0] = nvt[0];
                    evid[af][3-ii-jj] = nvt[0];
                }
            }
        }

        if(bf>=0 && lac>edge_thr)
        {
            int ii = faces[bf].indexof(a);
            int jj = faces[bf].indexof(c);
            int ov = faces[bf][3-ii-jj];
            if(bflag[ov])
            {
                int id = evid[f][1];
                if(id>=0)
                    nvt[1] = id;
                else
                {
                   // verts.push_back((va+vc)*0.5f);
                    verts.push_back((va+vc)*0.375f+(verts[b]+verts[ov])*0.125f);
                    bflag.push_back(true); bvertflag.push_back(true);
                    nvt[1] = verts.size()-1;
                    evid[f][1] = nvt[1];
                    evid[bf][3-ii-jj] = nvt[1];
                }
            }
        }

        if(cf>=0 && lab>edge_thr)
        {
            int ii = faces[cf].indexof(a);
            int jj = faces[cf].indexof(b);
            int ov = faces[cf][3-ii-jj];
            if(bflag[ov])
            {
                int id = evid[f][2];
                if(id>=0)
                    nvt[2] = id;
                else
                {
                  //  verts.push_back((va+vb)*0.5f);
                    verts.push_back((va+vb)*0.375f+(verts[c]+verts[ov])*0.125f);
                    bflag.push_back(true); bvertflag.push_back(true);
                    nvt[2] = verts.size()-1;
                    evid[f][2] = nvt[2];
                    evid[cf][3-ii-jj] = nvt[2];
                }
            }
        }

        if(nvt[0]<0&&nvt[1]<0&&nvt[2]<0)
            continue;
        if(nvt[0]>=0&&nvt[1]>=0&&nvt[2]>=0)
        {
            tfaces[f][0] = a; tfaces[f][1] = nvt[2]; tfaces[f][2] = nvt[1];
            tfaces.push_back(Face(nvt[2], b, nvt[0]));
            tfaces.push_back(Face(nvt[0], c, nvt[1]));
            tfaces.push_back(Face(nvt[2], nvt[0], nvt[1]));
        }

        if(nvt[0]>=0&&nvt[1]>=0&&nvt[2]<0)
        {
            tfaces[f][0] = a; tfaces[f][1] = b; tfaces[f][2] = nvt[0];
            tfaces.push_back(Face(a, nvt[0], nvt[1]));
            tfaces.push_back(Face(nvt[0], c, nvt[1]));
        }

        if(nvt[0]<0&&nvt[1]>=0&&nvt[2]>=0)
        {
            tfaces[f][0] = a; tfaces[f][1] = nvt[2]; tfaces[f][2] = nvt[1];
            tfaces.push_back(Face(nvt[1], nvt[2], c));
            tfaces.push_back(Face(nvt[2], b, c));
        }

        if(nvt[0]>=0&&nvt[1]<0&&nvt[2]>=0)
        {
            tfaces[f][0] = a; tfaces[f][1] = nvt[0]; tfaces[f][2] = c;
            tfaces.push_back(Face(a, nvt[2], nvt[0]));
            tfaces.push_back(Face(nvt[2], b, nvt[0]));
        }

        if(nvt[0]<0&&nvt[1]<0&&nvt[2]>=0)
        {
            tfaces[f][0] = a; tfaces[f][1] = nvt[2]; tfaces[f][2] = c;
            tfaces.push_back(Face(nvt[2], b, c));
        }

        if(nvt[0]>=0&&nvt[1]<0&&nvt[2]<0)
        {
            tfaces[f][0] = a; tfaces[f][1] = b; tfaces[f][2] = nvt[0];
            tfaces.push_back(Face(a, nvt[0], c));
        }

        if(nvt[0]<0&&nvt[1]>=0&&nvt[2]<0)
        {
            tfaces[f][0] = a; tfaces[f][1] = b; tfaces[f][2] = nvt[1];
            tfaces.push_back(Face(nvt[1], b, c));
        }
    }

    vts = verts; faces = tfaces;
    cout<<vts.size()<<" "<<faces.size()<<endl;
    fcs.clear();
    for(int f=0; f<faces.size(); f++)
    {
        fcs.push_back(faces[f][0]+1);
        fcs.push_back(faces[f][1]+1);
        fcs.push_back(faces[f][2]+1);
    }
}

