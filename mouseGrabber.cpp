#include "mouseGrabber.h"
#include <QGLViewer/manipulatedFrame.h>
#include <QGLWidget>
#include <ctime>
#include <cmath>
#include <fstream>
#include <array>
#include <iostream>
#include <fstream>
#include <igl/upsample.h>
#include <igl/per_vertex_attribute_smoothing.h>
#include "folder.h"
#include <GL/glu.h>

#define SetLightParamsV(L,ambt,difu,spec, shns){glLightfv(L,GL_AMBIENT,ambt);glLightfv(L,GL_DIFFUSE,difu);glLightfv(L,GL_SPECULAR,spec);glLightf(L,GL_SHININESS,shns);}
#define SetMaterialParamsV(F,ambt,difu,spec){glMaterialfv(F,GL_AMBIENT,ambt);glMaterialfv(F,GL_DIFFUSE,difu);glMaterialfv(F,GL_SPECULAR,spec);}
#define GetMaterialParamsV(F,ambt,difu,spec){glGetMaterialfv(F,GL_AMBIENT,ambt);glGetMaterialfv(F,GL_DIFFUSE,difu);glGetMaterialfv(F,GL_SPECULAR,spec);}

using namespace qglviewer;
using namespace std;

inline float distPoint(int x0, int y0, int x1, int y1)
{
    return sqrtf((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
}

inline void Set4Array(GLfloat* array, float a, float b, float c, float d){array[0] = a; array[1] = b; array[2] = c; array[3] = d;}

bool Viewer::loadfixedId(const char *path)
{
   fixedid.clear();
   ifstream fid(path);
   for(int i=0; i<264; i++)
   {
       int id; fid>>id;
       fixedid.push_back(id);
   }
   fid.close();
   return true;
}

bool Viewer::loadsidefixedId(const char* path)
{
    sidefixid.clear();
    ifstream fid(path);
    for(int i=0; i<1434; i++)
    {
        int id; fid>>id;
        sidefixid.push_back(id);
    }
    fid.close();
    return true;
}

bool Viewer::loadOBJ(const char * path){

        FILE * file = fopen(path, "r");
        if( file == NULL ){
                printf("Impossible to open the file ! Are you in the right path ? See Tutorial 1 for details\n");
                return false;
        }

        vector<int> new_vertexIndices;
        vts_label.resize(VERTEX_SIZE);

        while(1){
                char lineHeader[256];
                // read the first word of the line
                int res = fscanf(file, "%s", lineHeader);
                if (res == EOF){
                    printf("End Reading.\n");
                    break;
                }
//                if (!strcmp(lineHeader, "v")){
//                    float f1, f2, f3;
//                    fscanf(file, "%f %f %f\n", &f1, &f2, &f3);
//                    new_raw_vertices.push_back(QVector3D(f1, f2, f3));
//                }
                else if (!strcmp(lineHeader, "f")){
                    // case: f 1 2 3
                    int f1, f2, f3;
                    int t1, t2, t3;
                    int n1, n2, n3;
                    fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &f1, &t1, &n1, &f2, &t2, &n2, &f3, &t3, &n3);
                    new_vertexIndices.push_back(f1);
                    new_vertexIndices.push_back(f2);
                    new_vertexIndices.push_back(f3);
                    vts_label[f1-1] = t1-1;
                    vts_label[f2-1] = t2-1;
                    vts_label[f3-1] = t3-1;
                }  else if (!strcmp(lineHeader, "vt")){
                    // case: vt 1 2
                    float t1, t2;
                    fscanf(file, "%f %f\n", &t1, &t2);
                    vts.push_back(make_pair(t1, t2));
                } else {
                    // Probably a comment, eat up the rest of the line
                    char stupidBuffer[1000];
                    fgets(stupidBuffer, 1000, file);
                }
        }

        loadfixedId(std::string(localfolder()+"\\models\\\indices\\idfixed.in").c_str());
        loadsidefixedId(std::string(localfolder()+"\\models\\\indices\\idsidefixed.in").c_str());
        init_vertexIndices = new_vertexIndices;
        return true;

}

void Viewer::centerize(vector<vec3>& verts)
{
    float minx, miny, minz, maxx, maxy, maxz;
    minx = miny = minz = 1e+8;
    maxx = maxy = maxz = -1e+8;
    for(int i=0; i<verts.size(); i++)
    {
        vec3 vt = verts[i];
        if(vt[0]>maxx)  maxx = vt[0];
        if(vt[0]<minx)  minx = vt[0];
        if(vt[1]>maxy)  maxy = vt[1];
        if(vt[1]<miny)  miny = vt[1];
        if(vt[2]>maxz)  maxz = vt[2];
        if(vt[2]<minz)  minz = vt[2];
    }

    float x0, y0, z0, s;
    x0 = (minx+maxx)/2;
    y0 = (miny+maxy)/2;
    z0 = (minz+maxz)/2;

    s = 2.0/(maxx-minx);

    for(int i=0; i<verts.size(); i++)
    {
        vec3 vt = verts[i];
        vt[0] = (vt[0]-x0)*s;
        vt[1] = (vt[1]-y0)*s;
        vt[2] = (vt[2]-z0)*s;
        verts[i] = vt;
    }

}

bool Viewer::reloadOBJ(const double* new_vertices){

    ++currentpos;
    float current_max_dist = 0;
    for(int i=0; i< VERTEX_SIZE; ++i){
        if (backpos.find(i) != backpos.end()){
            float current_dist = new_vertices[3*i]*new_vertices[3*i] + new_vertices[3*i+1]*new_vertices[3*i+1] + new_vertices[3*i+2]*new_vertices[3*i+2];
            current_max_dist = max(current_max_dist, current_dist);
        }
    }
    current_scale = abs(max_dist/sqrt(current_max_dist));

    vector<QVector3D> new_raw_vertices;
    new_raw_vertices.resize(VERTEX_SIZE);
    // reset raw vetices
    for (int i=0; i < VERTEX_SIZE; ++i){
        new_raw_vertices[i] = QVector3D(new_vertices[3*i]*current_scale, new_vertices[3*i+1]*current_scale, new_vertices[3*i+2]*current_scale);
    }

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(init_vertexIndices);
    calcNeighbors();

    return reloadVertices();
}

bool Viewer::upsample(){

    int nvert = raw_vertices[currentpos].size();
    int nface = init_vertexIndices.size()/3;

    Eigen::MatrixXd v, new_v;
    Eigen::MatrixXi f;
    v.resize(nvert, 3); f.resize(nface, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = raw_vertices[currentpos][i].x();
        v(i, 1) = raw_vertices[currentpos][i].y();
        v(i, 2) = raw_vertices[currentpos][i].z();
    }

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = init_vertexIndices[3*i+j] - 1;
        }
    }

    igl::upsample(v, f, 1);
    igl::per_vertex_attribute_smoothing(v, f, new_v);

    vector<QVector3D> new_raw_vertices;
    vector<int> new_vertexIndices;
    int new_nvert = new_v.rows();
    int new_nface = f.rows();

    for (int i = 0; i < new_nvert; ++i){
        new_raw_vertices.push_back(QVector3D(new_v(i, 0), new_v(i, 1), new_v(i, 2)));
    }

    for (int i = 0; i < new_nface; ++i){
        new_vertexIndices.push_back(f(i, 0)+1);
        new_vertexIndices.push_back(f(i, 1)+1);
        new_vertexIndices.push_back(f(i, 2)+1);
    }

    raw_vertices[currentpos] = new_raw_vertices;
    //vertexIndices.push_back(new_vertexIndices);
    vertexIndices[currentpos] = new_vertexIndices;

    return true;
}

float Viewer::distThreshold()
{
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    GLdouble x0, x1, x2, x3, y0, y1, y2, y3, z;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    gluProject(1, -1, 1, modelview.data(), projection.data(), viewport.data(),
        &x0, &y0, &z);
    gluProject(-1, -1, 1, modelview.data(), projection.data(), viewport.data(),
        &x1, &y1, &z);
    gluProject(1, 1, 1, modelview.data(), projection.data(), viewport.data(),
        &x2, &y2, &z);
    gluProject(1, -1, -1, modelview.data(), projection.data(), viewport.data(),
        &x3, &y3, &z);
    vec2 a, b, c, d;
    a[0] = x0; a[1] = y0; b[0] = x1; b[1] = y1;
    c[0] = x2; c[1] = y2; d[0] = x3; d[1] = y3;

    float lab, lac, lad, len;
    lab = (a-b).length(); lac = (a-c).length(); lad = (a-d).length();
    len = lab>lac?lab:lac;
    if(lad>len) len = lad;
    return len;
}

void Viewer::smoothSampleArea()
{
      for(int iter= 0; iter<1; iter++)
      {
          // smoothing mesh
          for(int i=0; i<raw_vertices[currentpos].size(); i++)
          {
              if(bflag[i])
              {
                  vec3 vt(0, 0, 0); int ncount = 0;
                  for(std::set<int>::iterator it = neighbors[currentpos][i].begin();
                      it!=neighbors[currentpos][i].end(); ++it)
                  {
                      int jd = *it;
                      vec3 v = raw_vertices[currentpos][jd];
                      for(int k=0; k<3; k++)
                          vt[k]+=v[k];
                      ncount++;
                  }
                  if(ncount>0)
                  {
                      for(int k=0; k<3; k++)
                           raw_vertices[currentpos][i][k]=vt[k]/ncount;
                  }
              }
          }
      }

      reloadVertices();
}

void Viewer::smoothSampleAreaNormal()
{
      for(int iter= 0; iter<1; iter++)
      {
          // smoothing mesh
          for(int i=0; i<raw_vertices[currentpos].size(); i++)
          {
              if(bflag[i])
              {
                  vec3 nm(0, 0, 0); int ncount = 0;
                  for(std::set<int>::iterator it = neighbors[currentpos][i].begin();
                      it!=neighbors[currentpos][i].end(); ++it)
                  {
                      int jd = *it;
                      vec3 normal = raw_point_normals[currentpos][jd];
                      for(int k=0; k<3; k++)
                          nm[k]+=normal[k];
                      ncount++;
                  }
                  if(ncount>0)
                  {
                      for(int k=0; k<3; k++)
                           raw_point_normals[currentpos][i][k]=nm[k]/ncount;
                  }
              }
          }
      }
}

void Viewer::upsampleArea(int left, int right, int up, int down, vector<QPoint>& line){
    if (sampled){
        --currentpos;
        removeAfterCurrentpos();
    }

    int np = line.size();
    for (int i = 0; i + 1 < np; ++i){
        QPoint a = line[i];
        QPoint b = line[i+1];
        // Too Narrow
        float draw_dis = (a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y());
        if (draw_dis < 10 && i + 2 != np){
            line.erase(line.begin()+i+1, line.begin()+i+2);
            --i;
            --np;
        }
        // Too Wide
        else if (draw_dis > 50){
            line.insert(line.begin()+i+1, QPoint((a.x() + b.x())/2, (a.y() + b.y())/2));
            --i;
            ++np;
        }
    }
    int nq = raw_vertices[currentpos].size();
    int n_face = vertexIndices[currentpos].size();

    int width_extend = (right-left)/4 < 50 ? 50 : (right-left)/4;
    int height_extend = (down-up)/4 < 50 ? 50 : (down-up)/4;
    left < width_extend? left = 0: left -=width_extend;
    up < height_extend? up = 0: up -=height_extend;
    right > VIEW_SIZE - width_extend? right = VIEW_SIZE: right += width_extend;
    down > VIEW_SIZE - height_extend? down = VIEW_SIZE: down += height_extend;

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    const QVector3D viewdir (camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z);

    bflag.clear();
    bflag.resize(nq, false);

    float dist_thr = distThreshold()*0.05;

    cout<<"hanhan "<<dist_thr<<endl;

    for (int i = 0; i < nq; ++i){
        float cur_angle = QVector3D::dotProduct(viewdir, raw_point_normals[currentpos][i]);
        if(cur_angle < 0){ // the front points
            gluProject(raw_vertices[currentpos][i].x(), raw_vertices[currentpos][i].y(), raw_vertices[currentpos][i].z(),
                modelview.data(), projection.data(), viewport.data(),
                &current_x, &current_y, &current_z);
            current_y = VIEW_SIZE - current_y;
            int int_current_x = current_x, int_current_y = current_y;
            if (int_current_x >= left&& int_current_x <= right && int_current_y >= up && int_current_y <= down) {
                float dist = 1e+5;
                for(int k=0; k<line.size(); k++)
                {
                    float dis = distPoint(line[k].x(), line[k].y(), int_current_x, VIEW_SIZE-int_current_y);
                    if(dis<dist) dist = dis;
                }
                if(dist<dist_thr)
                    bflag[i] = true;
            }
        }
    }

    vector<int> pre_vertexIndices = vertexIndices[currentpos];
    vector<QVector3D> new_raw_vertices = raw_vertices[currentpos];
    vector<int> new_vertexIndices;

    map<pair<int, int>, int> new_vertex_map; // map <line start index, line end index> to new mid index
    vector<int> two_sample_face;
    vector<int> three_sample_face;
    vector<int> zero_sample_face;

    for(int it=0; it<2; it++)
    {
        new_vertexIndices.clear();
        n_face = pre_vertexIndices.size();


        // determine if a face is dense enough
        for(int i=0; i<n_face; i+=3)
        {
            int va = pre_vertexIndices[i]-1, vb = pre_vertexIndices[i+1]-1, vc = pre_vertexIndices[i+2]-1;
            if(!bflag[va]&&!bflag[vb]&&!bflag[vc]) continue;

            vec3 a = new_raw_vertices[va];
            vec3 b = new_raw_vertices[vb];
            vec3 c = new_raw_vertices[vc];

            if((a-b).length()<1e-2&&(a-c).length()<1e-2&&(c-b).length()<1e-2)
            {
                bflag[va] = bflag[vb] = bflag[vc] = false;
            }
        }




        new_vertex_map.clear();
        two_sample_face.clear();
        three_sample_face.clear();
        zero_sample_face.clear();

        for (int i = 0; i < n_face; i+= 3){
            int va = pre_vertexIndices[i]-1, vb = pre_vertexIndices[i+1]-1, vc = pre_vertexIndices[i+2]-1;

            if (bflag[va]&&bflag[vb]&&bflag[vc]){
                three_sample_face.push_back(va); three_sample_face.push_back(vb); three_sample_face.push_back(vc);
            } else if (!bflag[va]&&bflag[vb]&&bflag[vc])
            {
                two_sample_face.push_back(vb); two_sample_face.push_back(vc); two_sample_face.push_back(va);
            } else if(bflag[va]&&!bflag[vb]&&bflag[vc])
            {
                two_sample_face.push_back(vc); two_sample_face.push_back(va); two_sample_face.push_back(vb);
            } else if(bflag[va]&&bflag[vb]&&!bflag[vc])
            {
                 two_sample_face.push_back(va); two_sample_face.push_back(vb); two_sample_face.push_back(vc);
            }
            else {
                zero_sample_face.push_back(va); zero_sample_face.push_back(vb); zero_sample_face.push_back(vc);
            }
        }


        int line_startindex, line_endindex, first_mid_index, second_mid_index, third_mid_index;
        int first_quarter_index, first_midquarter_index, second_quarter_index, second_midquarter_index, third_quarter_index, third_midquarter_index,
                middle_first_index, middle_second_index, middle_third_index;
        pair<int, int> start_end, end_start;

        // for three_sample_face
        int nthree = three_sample_face.size();
        for (int i = 0; i < nthree; i+=3){

            // 1-level subdivision
            int first_index = three_sample_face[i], second_index = three_sample_face[i+1], third_index = three_sample_face[i+2];
            // first middle point
            line_startindex = first_index;
            line_endindex = second_index;
            start_end = make_pair(line_startindex, line_endindex);
            end_start = make_pair(line_endindex, line_startindex);
            if (new_vertex_map.find(start_end) == new_vertex_map.end()){
                new_raw_vertices.push_back(QVector3D(
                        (new_raw_vertices[line_startindex].x()+ new_raw_vertices[line_endindex].x())/2,
                        (new_raw_vertices[line_startindex].y()+ new_raw_vertices[line_endindex].y())/2,
                        (new_raw_vertices[line_startindex].z()+ new_raw_vertices[line_endindex].z())/2));
                bflag.push_back(true);
                first_mid_index = new_vertex_map[end_start] = new_vertex_map[start_end] = new_raw_vertices.size();
            } else first_mid_index = new_vertex_map[end_start];

            // second middle point
            line_startindex = second_index;
            line_endindex = third_index;
            start_end = make_pair(line_startindex, line_endindex);
            end_start = make_pair(line_endindex, line_startindex);
            if (new_vertex_map.find(start_end) == new_vertex_map.end()){
                new_raw_vertices.push_back(QVector3D(
                        (new_raw_vertices[line_startindex].x()+ new_raw_vertices[line_endindex].x())/2,
                        (new_raw_vertices[line_startindex].y()+ new_raw_vertices[line_endindex].y())/2,
                        (new_raw_vertices[line_startindex].z()+ new_raw_vertices[line_endindex].z())/2));
                bflag.push_back(true);
                second_mid_index = new_vertex_map[end_start] = new_vertex_map[start_end] = new_raw_vertices.size();
            } else second_mid_index = new_vertex_map[end_start];

            // third middle point
            line_startindex = third_index;
            line_endindex = first_index;
            start_end = make_pair(line_startindex, line_endindex);
            end_start = make_pair(line_endindex, line_startindex);
            if (new_vertex_map.find(start_end) == new_vertex_map.end()){
                new_raw_vertices.push_back(QVector3D(
                        (new_raw_vertices[line_startindex].x()+ new_raw_vertices[line_endindex].x())/2,
                        (new_raw_vertices[line_startindex].y()+ new_raw_vertices[line_endindex].y())/2,
                        (new_raw_vertices[line_startindex].z()+ new_raw_vertices[line_endindex].z())/2));
                bflag.push_back(true);

                third_mid_index = new_vertex_map[end_start] = new_vertex_map[start_end] = new_raw_vertices.size();
            } else third_mid_index = new_vertex_map[end_start];

            new_vertexIndices.push_back(first_index+1);
            new_vertexIndices.push_back(first_mid_index);
            new_vertexIndices.push_back(third_mid_index);

            new_vertexIndices.push_back(second_index+1);
            new_vertexIndices.push_back(second_mid_index);
            new_vertexIndices.push_back(first_mid_index);

            new_vertexIndices.push_back(third_index+1);
            new_vertexIndices.push_back(third_mid_index);
            new_vertexIndices.push_back(second_mid_index);

            new_vertexIndices.push_back(first_mid_index);
            new_vertexIndices.push_back(second_mid_index);
            new_vertexIndices.push_back(third_mid_index);

         }

        // for two_sample_face
        int ntwo = two_sample_face.size();
        for (int i = 0; i < ntwo; i+=3){
            int first_index = two_sample_face[i], second_index = two_sample_face[i+1], third_index = two_sample_face[i+2];
            // first middle point
            line_startindex = first_index;
            line_endindex = second_index;
            start_end = make_pair(line_startindex, line_endindex);
            end_start = make_pair(line_endindex, line_startindex);
            if (new_vertex_map.find(start_end) != new_vertex_map.end()){
                first_mid_index = new_vertex_map[end_start];
                new_vertexIndices.push_back(first_index+1);
                new_vertexIndices.push_back(first_mid_index);
                new_vertexIndices.push_back(third_index+1);

                new_vertexIndices.push_back(first_mid_index);
                new_vertexIndices.push_back(second_index+1);
                new_vertexIndices.push_back(third_index+1);
            } else {
                new_vertexIndices.push_back(first_index+1);
                new_vertexIndices.push_back(second_index+1);
                new_vertexIndices.push_back(third_index+1);
            }
        }

        // for zero_sample_face
        int nzero = zero_sample_face.size();
        for (int i = 0; i < nzero; ++i){
            new_vertexIndices.push_back(zero_sample_face[i]+1);
        }

        pre_vertexIndices = new_vertexIndices;
    }

    ++currentpos;
    vertexIndices.push_back(new_vertexIndices);
    raw_vertices.push_back(new_raw_vertices);
    for (int i = hightlight_deformneighbor_rawindex.size(); i < raw_vertices[currentpos].size(); ++i){
        hightlight_deformneighbor_rawindex.push_back(false);
        hightlight_selection_rawindex.push_back(false);
    }

    calcNeighbors();

    smoothSampleArea();

    //reloadVertices();

    update();
    sampled = true;

}

void Viewer::upsamplelinearea(int left, int right, int up, int down, vector<QPoint> &line)
{
    if (sampled){
        --currentpos;
        removeAfterCurrentpos();
    }

    int np = line.size();
    for (int i = 0; i + 1 < np; ++i){
        QPoint a = line[i];
        QPoint b = line[i+1];
        // Too Narrow
        float draw_dis = (a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y());
        if (draw_dis < 8 && i + 2 != np){
            line.erase(line.begin()+i+1, line.begin()+i+2);
            --i;
            --np;
        }
        // Too Wide
        else if (draw_dis > 40){
            line.insert(line.begin()+i+1, QPoint((a.x() + b.x())/2, (a.y() + b.y())/2));
            --i;
            ++np;
        }
    }

    int nq = raw_vertices[currentpos].size();
    int n_face = vertexIndices[currentpos].size();

    int width_extend = (right-left)/4 < 50 ? 50 : (right-left)/4;
    int height_extend = (down-up)/4 < 50 ? 50 : (down-up)/4;
    left < width_extend? left = 0: left -=width_extend;
    up < height_extend? up = 0: up -=height_extend;
    right > VIEW_SIZE - width_extend? right = VIEW_SIZE: right += width_extend;
    down > VIEW_SIZE - height_extend? down = VIEW_SIZE: down += height_extend;

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    const QVector3D viewdir (camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z);

    vec3 campos(camera()->position().x, camera()->position().y, camera()->position().z);


    bflag.clear(); bflag.resize(nq, false);
    float dist_thr = 0;

    vector< vector<float>> dists;
    dists.resize(line.size());

    for(int k=0; k<line.size(); k++)
    {
        dists[k].resize(nq, -1);
        float dist = 1e+5; int id = -1;
        for (int i = 0; i < nq; ++i){
            float cur_angle = QVector3D::dotProduct(viewdir, raw_point_normals[currentpos][i]);
            if(cur_angle < 0){ // the front points
                gluProject(raw_vertices[currentpos][i].x(), raw_vertices[currentpos][i].y(), raw_vertices[currentpos][i].z(),
                    modelview.data(), projection.data(), viewport.data(),
                    &current_x, &current_y, &current_z);
                current_y = VIEW_SIZE - current_y;
                int int_current_x = current_x, int_current_y = current_y;
                if (int_current_x >= left&& int_current_x <= right && int_current_y >= up && int_current_y <= down) {
                     float dis = distPoint(line[k].x(), line[k].y(), int_current_x, VIEW_SIZE-int_current_y);
                     dists[k][i] = dis;
                     if(dis<dist)
                     {
                        dist = dis; id = i;
                     }
                }
            }
        }

        if(id>=0)
        {
            float len = 0; int ncount = 0;
            vec3 p = raw_vertices[currentpos][id];
            GLdouble x0, y0, z0, x1, y1, z1;
            gluProject(p.x(), p.y(), p.z(),modelview.data(), projection.data(), viewport.data(), &x0, &y0, &z0);

            for(auto it = neighbors[currentpos][id].begin(); it!=neighbors[currentpos][id].end(); it++)
            {
               vec3 q = raw_vertices[currentpos][*it];
               gluProject(q.x(), q.y(), q.z(),modelview.data(), projection.data(), viewport.data(), &x1, &y1, &z1);
               len+=sqrtf((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)); ncount++;
            }

            if(ncount>0) len/=ncount;
            if(len>dist_thr) dist_thr = len ;
        }
    }

    for(int i=0; i<nq; i++)
    {
        float dist = 1e+5;
        for(int k=0; k<line.size(); k++)
        {
            if(dists[k][i]>=0&&dists[k][i]<dist)
                dist = dists[k][i];
        }
        if(dist<dist_thr*3)
            bflag[i] = true;
    }

    vector<bool> bsampleflag = bflag;

    vector<vec3> new_raw_vertices = raw_vertices[currentpos];
    vector<int> new_vertexIndices = vertexIndices[currentpos];

    cout<<"start upsampe"<<endl;
    subd->upsample(new_raw_vertices, new_vertexIndices, bsampleflag, bflag, 1.5e-2);

    for(int i=0; i<bsampleflag.size(); i++)
    {
        if(bsampleflag[i])
        {
            float dist = 1e+5;
            for(int k=0; k<line.size(); k++)
            {
                gluProject(new_raw_vertices[i].x(), new_raw_vertices[i].y(), new_raw_vertices[i].z(),
                    modelview.data(), projection.data(), viewport.data(),
                    &current_x, &current_y, &current_z);
                int int_current_x = current_x, int_current_y = current_y;

                float dis = distPoint(line[k].x(), line[k].y(), int_current_x, int_current_y);
                if(dis<dist)
                 {
                     dist = dis;
                 }

            }
            if(dist>dist_thr*1.5)
                bsampleflag[i] = false;
        }
    }

    subd->upsample(new_raw_vertices, new_vertexIndices, bsampleflag, bflag, 1.5e-2);

    for(int i=0; i<bsampleflag.size(); i++)
    {
        if(bsampleflag[i])
        {
            float dist = 1e+5;
            for(int k=0; k<line.size(); k++)
            {
                gluProject(new_raw_vertices[i].x(), new_raw_vertices[i].y(), new_raw_vertices[i].z(),
                    modelview.data(), projection.data(), viewport.data(),
                    &current_x, &current_y, &current_z);
                int int_current_x = current_x, int_current_y = current_y;

                float dis = distPoint(line[k].x(), line[k].y(), int_current_x, int_current_y);
                if(dis<dist)
                 {
                     dist = dis;
                 }

            }
            if(dist>dist_thr*1)
                bsampleflag[i] = false;
        }
    }

    subd->upsample(new_raw_vertices, new_vertexIndices, bsampleflag, bflag, 1.5e-2);

    ++currentpos;
    vertexIndices.push_back(new_vertexIndices);
    raw_vertices.push_back(new_raw_vertices);


    calcNeighbors();

    reloadVertices();
    smoothSampleAreaNormal();
    update();
    sampled = true;
}

bool Viewer::reloadVertices(){
    int new_nvert = raw_vertices[currentpos].size();
    int new_nindex = vertexIndices[currentpos].size();

    // check upsample size change
    vector<QVector3D> new_vertices; new_vertices.resize(new_nindex);
    vector<QVector3D> new_poly_normals; new_poly_normals.resize(new_nindex/3);
    vector<QVector3D> new_point_normals; new_point_normals.resize(new_nindex);
    vector<QVector3D> new_raw_point_normals; new_raw_point_normals.resize(new_nvert);

    // start reloading
    for(int i=0; i<new_nindex; i++){
        new_vertices[i] = raw_vertices[currentpos][vertexIndices[currentpos][i]-1];
    }
    vertices.push_back(new_vertices);

    vector<vector<QVector3D> > vertex_maps_normals;
    vertex_maps_normals.resize(new_nvert);

    for (int i=0; i<new_nindex; i=i+3){
        QVector3D Vu = new_vertices[i+1] - new_vertices[i];
        QVector3D Vv = new_vertices[i+2] - new_vertices[i];
        QVector3D surfaceNormal(Vu.y() * Vv.z() - Vu.z() * Vv.y(),
                               Vu.z() * Vv.x() - Vu.x() * Vv.z(),
                               Vu.x() * Vv.y() - Vu.y() * Vv.x());
        surfaceNormal.normalize();

        new_poly_normals[i/3] = surfaceNormal;

        vertex_maps_normals[vertexIndices[currentpos][i]-1].push_back(surfaceNormal);
        vertex_maps_normals[vertexIndices[currentpos][i+1]-1].push_back(surfaceNormal);
        vertex_maps_normals[vertexIndices[currentpos][i+2]-1].push_back(surfaceNormal);
    }

    poly_normals.push_back(new_poly_normals);

    // take mean normal of vertices
    for (int i=0; i<new_nvert;++i){
        QVector3D basetemp(0, 0, 0);
        for (int j=0; j<vertex_maps_normals[i].size();++j){
            basetemp += vertex_maps_normals[i][j];
        }
        basetemp.normalize();
        new_raw_point_normals[i] = basetemp;
    }
    raw_point_normals.push_back(new_raw_point_normals);

    // map vertex normals to index normals
    for(int i=0; i<new_nindex; i++){
        new_point_normals[i] = new_raw_point_normals[vertexIndices[currentpos][i]-1];
    }
    point_normals.push_back(new_point_normals);
    resizeVectors();

    return true;
}

void Viewer::smoothlinevertices(vector<int>& line)
{
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<GLdouble> lx, ly, lz;
    lx.resize(line.size(), 0);
    ly.resize(line.size(), 0);
    lz.resize(line.size(), 0);

    for(int i=0; i<line.size(); i++)
    {
        int id = line[i];
        std::array<GLdouble, 3> screen_coords;
        gluProject(raw_vertices[currentpos][id].x(), raw_vertices[currentpos][id].y(), raw_vertices[currentpos][id].z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        lx[i] = screen_coords[0]; ly[i] = screen_coords[1]; lz[i] = screen_coords[2];
    }

    // smoothing
    for(int it=0; it<10; it++)
    {
       for(int i=1; i<lx.size()-1; i++)
       {
           lx[i]= (lx[i-1]+lx[i+1])*0.5f;
           ly[i]= (ly[i-1]+ly[i+1])*0.5f;
       }
    }

    for(int i= 0; i<lx.size(); i++)
    {
        GLdouble newobj_vertex[3];
        gluUnProject(lx[i], ly[i], lz[i], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);
        int id = line[i];
        raw_vertices[currentpos][id][0] = newobj_vertex[0];
        raw_vertices[currentpos][id][1] = newobj_vertex[1];
        raw_vertices[currentpos][id][2] = newobj_vertex[2];
    }

}

void Viewer::calshiftvts(vector<int>& line, vector<vec3>& vts,  float zdir)
{
    int np = line.size();

    cout<<"han han "<<np<<endl;

    vector<float> wt;
    wt.resize(np, 1.0);
    wt[0] = wt[np-1] = 0.0f;
    wt[1] = wt[np-2] = 0.5f;

    for(int it = 0; it<5; it++)
    {
        for(int i=1; i<np-1; i++)
        {
            wt[i] = (wt[i-1]+wt[i+1])*0.5f;
        }
    }

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<GLdouble> lx, ly, lz;
    lz.resize(np, 0); lx.resize(np, 0); ly.resize(np, 0);

    for (int i = 0; i <np; ++i){
        float s = wt[i];
        int id = line[i];

        std::array<GLdouble, 3> screen_coords;
        gluProject(raw_vertices[currentpos][id].x(), raw_vertices[currentpos][id].y(), raw_vertices[currentpos][id].z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        lx[i] = screen_coords[0]; ly[i] = screen_coords[1]; lz[i] = screen_coords[2]-s*zdir;
    }

    for(int i=0; i<np; i++)
    {
        GLdouble newobj_vertex[3];
        gluUnProject(lx[i], ly[i], lz[i], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        vts.push_back(QVector3D(newobj_vertex[0],newobj_vertex[1], newobj_vertex[2]));
    }


}

void Viewer::resizeVectors(){
    int to_erase = currentpos > 15 ? currentpos - 15 : 0;
    currentpos -= to_erase;
    vertexIndices.erase(vertexIndices.begin(), vertexIndices.begin()+to_erase);
    raw_vertices.erase(raw_vertices.begin(), raw_vertices.begin()+to_erase);
    vertices.erase(vertices.begin(), vertices.begin()+to_erase);
    raw_point_normals.erase(raw_point_normals.begin(), raw_point_normals.begin()+to_erase);
    point_normals.erase(point_normals.begin(), point_normals.begin()+to_erase);
    poly_normals.erase(poly_normals.begin(), poly_normals.begin()+to_erase);
    neighbors.erase(neighbors.begin(), neighbors.begin()+to_erase);
}

void Viewer::init(){

    // load texture
    textureImg = QGLWidget::convertToGLFormat(QImage(QString::fromStdString(localfolder()+"\\models\\texture\\texture.png")));
    glTexImage2D(GL_TEXTURE_2D, 0, 4, textureImg.width(), textureImg.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, textureImg.bits());

    // Absolutely needed for MouseGrabber
    setMouseTracking(true);

    // In order to make the manipulatedFrame displacements clearer
    //setAxisIsDrawn();

    // Initialize the CameraPathPlayer MouseGrabber array
    nbPlayers_ = 12;
    player_ = new CameraPathPlayer*[nbPlayers_];
    for (int i=0; i<nbPlayers_; ++i)
        player_[i] = NULL;

    //load back vectors of the model - to set the bounding box
    ifstream fs(localfolder()+"\\models\\boundingBoxes\\backIndices.txt");
    int* backs = (int *)malloc(backnum*sizeof(int));
    for (int i = 0; i < backnum; ++i) fs >> backs[i];
    backpos.insert(backs, backs+backnum);
    free(backs);


    loadOBJ(std::string(localfolder()+"\\objs\\shape.obj").c_str());
    lapdeform = new LapDeform;
    subd = new sketchRender;
    restoreStateFromFile();
    setEnvironment(255, 255, 255); // three entries correspond to background color;

    int nvert = VERTEX_SIZE;

    for (int i = 0; i < nvert; ++i) hightlight_selection_rawindex.push_back(false);
    for (int i = 0; i < nvert; ++i) hightlight_deformneighbor_rawindex.push_back(false);
    for (int i = 0; i < nvert; ++i) hightlight_init.push_back(false);

    printf("viewer finished\n");
}

void Viewer::calc2Dshape(){

    // get screen info for 2D mapping
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    // reset shape vertices
    shape_vertices_2D.erase(shape_vertices_2D.begin(), shape_vertices_2D.end());
    shape_vertices_3D.erase(shape_vertices_3D.begin(), shape_vertices_3D.end());
    shape_vertices_index.erase(shape_vertices_index.begin(), shape_vertices_index.end());
    shape_vertices_2D_z_value.erase(shape_vertices_2D_z_value.begin(), shape_vertices_2D_z_value.end());

    const float cos_min_degree = cos(MIN_NORMAL_DEGREE* PI / 180.0);
    const float cos_max_degree = cos(MAX_NORMAL_DEGREE* PI / 180.0);

    // calculate vertices with normal within MIN_NORMAL_DEGREE ~ MAX_NORMAL_DEGREE
    const QVector3D viewdir (camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z);
    const float viewdir_len = viewdir.length();
    int nvert = raw_vertices[currentpos].size();
    for (int i = 0; i < nvert; ++i){
        QVector3D cur_normal = raw_point_normals[currentpos][i];
        float cur_angle = (QVector3D::dotProduct(viewdir, cur_normal))/(cur_normal.length()*viewdir_len);
        if(cur_angle > cos_max_degree && cur_angle < cos_min_degree){
            shape_vertices_3D.push_back(raw_vertices[currentpos][i]);
            shape_vertices_index.push_back(i);
            std::array<GLdouble, 3> screen_coords;
            gluProject(raw_vertices[currentpos][i].x(), raw_vertices[currentpos][i].y(), raw_vertices[currentpos][i].z(),
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);
            QVector2D cur_2D_vertex(screen_coords[0], VIEW_SIZE - screen_coords[1]);
            shape_vertices_2D_z_value.push_back(screen_coords[2]);
            shape_vertices_2D.push_back(cur_2D_vertex);
        }
    }

    printf("Total %d Vertices Selected.\n\n", shape_vertices_3D.size());

}

QVector3D Viewer::getCameraPos(){
    QVector3D pos;
    pos[0] = camera()->position().x;
    pos[1] = camera()->position().y;
    pos[2] = camera()->position().z;
    return pos;
}

void Viewer::projectfrontsketches(vector< vector<QVector3D>>& sketches, vector< vector<QPoint>>& lines)
{
    MoveToFront();
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    std::array<GLdouble, 3> screen_coords;

    // map to 2D
    for(int i=0; i<sketches.size(); i++)
    {
        vector<QPoint> line;
        for(int j=0; j<sketches[i].size(); j++)
        {
            gluProject(sketches[i][j][0], sketches[i][j][1], sketches[i][j][2],
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            QPoint p((int)screen_coords[0], (int)(VIEW_SIZE-screen_coords[1]));
            line.push_back(p);
        }
        // smoothing line
        {
            vector<QPoint> templine = line;
            for(int it=0; it<3; it++)
            {
                for(int i=1; i<line.size()-1; i++)
                {
                    QPoint a, b;
                    int x, y;
                    a = templine[i-1]; b = templine[i+1];
                    x = (int)((a.x()+b.x())/2.0f);
                    y = (int)((a.y()+b.y())/2.0f);
                    line[i] = QPoint(x, y);
                }
                templine = line;
            }
        }

        lines.push_back(line);
    }

    // reset camera
}

void Viewer::projectsidesketches(vector< vector<QVector3D>>& sketches, vector< vector<QPoint>>& lines)
{
    MoveToSide();
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    std::array<GLdouble, 3> screen_coords;

    // map to 2D
    for(int i=0; i<sketches.size(); i++)
    {
        vector<QPoint> line;
        for(int j=0; j<sketches[i].size(); j++)
        {
            gluProject(sketches[i][j][0], sketches[i][j][1], sketches[i][j][2],
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            QPoint p((int)screen_coords[0], (int)(VIEW_SIZE-screen_coords[1]));
            line.push_back(p);
        }
        // smoothing line
        {
            vector<QPoint> templine = line;
            for(int it=0; it<1; it++)
            {
                for(int i=1; i<line.size()-1; i++)
                {
                    QPoint a, b;
                    int x, y;
                    a = templine[i-1]; b = templine[i+1];
                    x = (int)((a.x()+b.x())/2.0f);
                    y = (int)((a.y()+b.y())/2.0f);
                    line[i] = QPoint(x, y);
                }
                templine = line;
            }
        }

        lines.push_back(line);
    }
}

void Viewer::printCurrentView()
{

}

void Viewer::projectLeftSideSketches(vector< vector<QVector3D>>& sketches, vector< vector<QPoint>>& lines)
{
    MoveToLeftSide();
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    std::array<GLdouble, 3> screen_coords;

    // map to 2D
    for(int i=0; i<sketches.size(); i++)
    {
        vector<QPoint> line;
        for(int j=0; j<sketches[i].size(); j++)
        {
            gluProject(sketches[i][j][0], sketches[i][j][1], sketches[i][j][2],
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            QPoint p((int)screen_coords[0], (int)(VIEW_SIZE-screen_coords[1]));
            line.push_back(p);
        }
        // smoothing line
        {
            vector<QPoint> templine = line;
            for(int it=0; it<3; it++)
            {
                for(int i=1; i<line.size()-1; i++)
                {
                    QPoint a, b;
                    int x, y;
                    a = templine[i-1]; b = templine[i+1];
                    x = (int)((a.x()+b.x())/2.0f);
                    y = (int)((a.y()+b.y())/2.0f);
                    line[i] = QPoint(x, y);
                }
                templine = line;
            }
        }

        lines.push_back(line);
    }
}

void Viewer::projectRightSideSketches(vector< vector<QVector3D>>& sketches, vector< vector<QPoint>>& lines)
{
    MoveToRightSide();
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    std::array<GLdouble, 3> screen_coords;

    // map to 2D
    for(int i=0; i<sketches.size(); i++)
    {
        vector<QPoint> line;
        for(int j=0; j<sketches[i].size(); j++)
        {
            gluProject(sketches[i][j][0], sketches[i][j][1], sketches[i][j][2],
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            QPoint p((int)screen_coords[0], (int)(VIEW_SIZE-screen_coords[1]));
            line.push_back(p);
        }
        // smoothing line
        {
            vector<QPoint> templine = line;
            for(int it=0; it<3; it++)
            {
                for(int i=1; i<line.size()-1; i++)
                {
                    QPoint a, b;
                    int x, y;
                    a = templine[i-1]; b = templine[i+1];
                    x = (int)((a.x()+b.x())/2.0f);
                    y = (int)((a.y()+b.y())/2.0f);
                    line[i] = QPoint(x, y);
                }
                templine = line;
            }
        }

        lines.push_back(line);
    }
}

void Viewer::calcNeighbors(){
    int nvert = raw_vertices[currentpos].size();
    int nface = vertexIndices[currentpos].size()/3;
    vector<set<int>> new_neighbors;
    new_neighbors.resize(nvert);
    for (int i = 0; i < nface; i++) {
        int va = vertexIndices[currentpos][3*i]-1;
        int vb = vertexIndices[currentpos][3*i+1]-1;
        int vc = vertexIndices[currentpos][3*i+2]-1;
        new_neighbors[va].insert(vb); new_neighbors[va].insert(vc);
        new_neighbors[vb].insert(va); new_neighbors[vb].insert(vc);
        new_neighbors[vc].insert(va); new_neighbors[vc].insert(vb);
    }
    neighbors.push_back(new_neighbors);
}

void Viewer::handleDeform(vector<pair<int, vec2>>& shiftmap, vector<bool>& bfixed)
{
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());


    vector<vec3> new_raw_vertices = raw_vertices[currentpos];

    vector<vec3> targetvt;
    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
    vector<int> count;
    count.resize(new_raw_vertices.size(), 0);

    for(int i=0; i<shiftmap.size(); i++)
    {
        int id = shiftmap[i].first;
        vec2 pt = shiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

        std::array<GLdouble, 3> screen_coords;
        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
        gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
        count[id]++;
    }

    vector<int> handleid;
    vector<vec3> handlevt;

    for(int i=0; i<fixedid.size(); i++)
    {
        int id = fixedid[i];
        handleid.push_back(id);
        handlevt.push_back(new_raw_vertices[id]);
        count[id] = -1;
    }

    for(int i=0; i<count.size(); i++)
    {
        if(count[i]<=0) continue;

        {
            handleid.push_back(i);
            vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);
            handlevt.push_back(tv);
        }
    }


    lapdeform->handle_deform(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
    //lapdeform->save_deform_info(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);


    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    --currentpos;
    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);

    ++currentpos;
    cout<<currentpos<<" "<<raw_vertices.size()<<endl;
    reloadVertices();
}

void Viewer::inittransform(vector<vec3>& verts, vector<int>& vid, vector<vec3>& targetvt)
{
    vec3 scent(0, 0, 0);
    vec3 tcent(0, 0, 0);
    for(int i=0; i<vid.size(); i++)
    {
        int id = vid[i];
        for(int k=0; k<3; k++)
        {
            scent[k]+=verts[id][k];
            tcent[k]+=targetvt[id][k];
        }
    }
    for(int k=0; k<3; k++)
    {
        scent[k]=scent[k]*(1.0f/vid.size());
        tcent[k]=tcent[k]*(1.0f/vid.size());
    }

    // cal scale and translation
    float len1 = 0;
    float len2 = 0;
    for(int i=0; i<vid.size(); i++)
    {
        int id = vid[i];
        vec3 os = verts[id]-scent;
        vec3 ot = targetvt[id]-tcent;
        len1 = len1+os.length();
        len2 = len2+ot.length();
    }
    len1 = len1*(1.0f/vid.size());
    len2 = len2*(1.0f/vid.size());
    float scale = len2/len1;
    // rescale
    for(int i=0; i<verts.size(); i++)
    {
        for(int j=0; j<3; j++)
        {
            verts[i][j] = scale*(verts[i][j]-scent[j])+ tcent[j];
        }
    }
}

void Viewer::handleDeform(vector<pair<int, vec2>>& shiftmap, const double* new_vertices, vector<bool>& bfixed)
{
    float current_max_dist = 0;
    for(int i=0; i< VERTEX_SIZE; ++i){
        if (backpos.find(i) != backpos.end()){
            float current_dist = new_vertices[3*i]*new_vertices[3*i] + new_vertices[3*i+1]*new_vertices[3*i+1] + new_vertices[3*i+2]*new_vertices[3*i+2];
            current_max_dist = max(current_max_dist, current_dist);
        }
    }
    current_scale = abs(max_dist/sqrt(current_max_dist));

    vector<vec3> new_raw_vertices;
    new_raw_vertices.resize(VERTEX_SIZE);
    // reset raw vetices
    for (int i=0; i < VERTEX_SIZE; ++i){
        new_raw_vertices[i] = vec3(new_vertices[3*i]*current_scale, new_vertices[3*i+1]*current_scale, new_vertices[3*i+2]*current_scale);
    }

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<vec3> cur_raw_vertices = raw_vertices[currentpos];
    vector<vec3> targetvt;
    targetvt.resize(cur_raw_vertices.size(), vec3(0, 0, 0));
    vector<int> count;
    count.resize(cur_raw_vertices.size(), 0);

    for(int i=0; i<shiftmap.size(); i++)
    {
        int id = shiftmap[i].first;
        vec2 pt = shiftmap[i].second;

       // cout<<id<<" "<<pt[0]<<" "<<pt[1]<<endl;
        float x, y, z;
        x = y = z = 0;

        std::array<GLdouble, 3> screen_coords;
        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
               modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
           gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                    projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        x = newobj_vertex[0], y = newobj_vertex[1], z = new_raw_vertices[id].z();
        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
        count[id]++;
    }

    vector<int> handleid;

    for(int i=0; i<count.size(); i++)
    {
      if(count[i]<=0) continue;
        handleid.push_back(i);
        vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);
       targetvt[i] = tv;
    }
    inittransform(new_raw_vertices, handleid, targetvt);


    targetvt.clear();
    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
    count.clear();
    count.resize(new_raw_vertices.size(), 0);

    for(int i=0; i<shiftmap.size(); i++)
    {
        int id = shiftmap[i].first;
        vec2 pt = shiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

            std::array<GLdouble, 3> screen_coords;
            gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            GLdouble newobj_vertex[3];
            gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                    projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

            x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
            targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
            count[id]++;

    }

    handleid.clear();
    vector<vec3> handlevt;

    for(int i=0; i<fixedid.size(); i++)
    {
        int id = fixedid[i];
        handleid.push_back(id);
        handlevt.push_back(cur_raw_vertices[id]);
        count[id] = -1;
    }

    for(int i=0; i<count.size(); i++)
    {
       if(count[i]<=0) continue;
       handleid.push_back(i);
            vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);

       handlevt.push_back(tv);
    }

    //lapdeform->save_deform_info(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
    lapdeform->handle_deform(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);

    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    --currentpos;
    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);

    ++currentpos;
    cout<<currentpos<<" "<<raw_vertices.size()<<endl;
    reloadVertices();
}

void Viewer::handleDeformAllView(vector<pair<int, vec2>>& frontshiftmap, vector<pair<int, vec2>>& sideshiftmap, const double* new_vertices,
                                 vector<bool>& bfrontfixed, vector<bool>& bsidefixed)
{
    float current_max_dist = 0;
    for(int i=0; i< VERTEX_SIZE; ++i){
        if (backpos.find(i) != backpos.end()){
            float current_dist = new_vertices[3*i]*new_vertices[3*i] + new_vertices[3*i+1]*new_vertices[3*i+1] + new_vertices[3*i+2]*new_vertices[3*i+2];
            current_max_dist = max(current_max_dist, current_dist);
        }
    }
    current_scale = abs(max_dist/sqrt(current_max_dist));

    vector<vec3> new_raw_vertices;
    new_raw_vertices.resize(VERTEX_SIZE);
    // reset raw vetices
    for (int i=0; i < VERTEX_SIZE; ++i){
        new_raw_vertices[i] = vec3(new_vertices[3*i]*current_scale, new_vertices[3*i+1]*current_scale, new_vertices[3*i+2]*current_scale);
    }

    MoveToFront();
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<vec3> cur_raw_vertices = raw_vertices[currentpos];
    vector<vec3> targetvt;
    targetvt.resize(cur_raw_vertices.size(), vec3(0, 0, 0));
    vector<int> count;
    count.resize(cur_raw_vertices.size(), 0);



    for(int i=0; i<frontshiftmap.size(); i++)
    {
        int id = frontshiftmap[i].first;
        vec2 pt = frontshiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

        std::array<GLdouble, 3> screen_coords;
        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
               modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
           gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                    projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        x = newobj_vertex[0], y = newobj_vertex[1], z = new_raw_vertices[id].z();
        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
        count[id]++;
    }

    vector<int> handleid;

    for(int i=0; i<count.size(); i++)
    {
      if(count[i]<=0) continue;
        handleid.push_back(i);
        vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);
       targetvt[i] = tv;
    }
    inittransform(new_raw_vertices, handleid, targetvt);


    targetvt.clear();
    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
    count.clear();
    count.resize(new_raw_vertices.size(), 0);

    vector<bool> bfixed; bfixed.resize(new_raw_vertices.size(), false);

    for(int i=0; i<bfrontfixed.size(); i++)
    {
        if(bfrontfixed[i]) bfixed[i] = true;
    }

    for(int i=0; i<bsidefixed.size(); i++)
    {
        if(bsidefixed[i]) bfixed[i] = true;
    }

    MoveToFront();
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    for(int i=0; i<frontshiftmap.size(); i++)
    {
        int id = frontshiftmap[i].first;
        vec2 pt = frontshiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

            std::array<GLdouble, 3> screen_coords;
            gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            GLdouble newobj_vertex[3];
            gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                    projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

            x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
            targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
            count[id]++;
    }

    MoveToSide();
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    for(int i=0; i<sideshiftmap.size(); i++)
    {
        int id = sideshiftmap[i].first;
        vec2 pt = sideshiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

            std::array<GLdouble, 3> screen_coords;
            gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
                modelview.data(), projection.data(), viewport.data(),
                screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

            GLdouble newobj_vertex[3];
            gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                    projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

            x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
            targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
            count[id]++;

    }

    handleid.clear();
    vector<vec3> handlevt;

    for(int i=0; i<fixedid.size(); i++)
    {
        int id = fixedid[i];
        handleid.push_back(id);
        handlevt.push_back(cur_raw_vertices[id]);
        count[id] = -1;
    }

    for(int i=0; i<count.size(); i++)
    {
       if(count[i]<=0) continue;
       handleid.push_back(i);
            vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);

       handlevt.push_back(tv);
    }

    //lapdeform->save_deform_info2(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
    lapdeform->handle_deform_all(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(vertexIndices[currentpos]);

    ++currentpos;
    reloadVertices();
}

//void Viewer::handleDeformSideView(vector<pair<int, vec2> > &shiftmap, vector<bool> &bfixed)
//{
//    std::array<GLdouble, 16> projection;
//    std::array<GLdouble, 16> modelview;
//    std::array<GLint, 4> viewport;
//    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
//    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
//    glGetIntegerv(GL_VIEWPORT, viewport.data());

//    vector<vec3> new_raw_vertices = raw_vertices[currentpos];

//    vector<vec3> targetvt;
//    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
//    vector<int> count;
//    count.resize(new_raw_vertices.size(), 0);

//    for(int i=0; i<shiftmap.size(); i++)
//    {
//        int id = shiftmap[i].first;
//        vec2 pt = shiftmap[i].second;

//        float x, y, z;
//        x = y = z = 0;

//        std::array<GLdouble, 3> screen_coords;
//        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
//            modelview.data(), projection.data(), viewport.data(),
//            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

//        GLdouble newobj_vertex[3];
//        gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
//                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

//        x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
//        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
//        count[id]++;
//    }

//    vector<int> handleid;
//    vector<vec3> handlevt;

//    for(int i=0; i<fixedid.size(); i++)
//    {
//        int id = fixedid[i];
//        handleid.push_back(id);
//        handlevt.push_back(new_raw_vertices[id]);
//        count[id] = -1;
//    }

//    for(int i=0; i<count.size(); i++)
//    {
//        if(count[i]<=0) continue;

//        {
//            handleid.push_back(i);
//            vec3 tv(0, 0, 0);
//            tv[0] = targetvt[i][0]/(1.0f*count[i]);
//             tv[1] = targetvt[i][1]/(1.0f*count[i]);
//              tv[2] = targetvt[i][2]/(1.0f*count[i]);
//            handlevt.push_back(tv);
//        }
//    }


//    //lapdeform->handle_deform(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
//    lapdeform->save_deform_info2(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
//    raw_vertices.push_back(new_raw_vertices);
//    vertexIndices.push_back(vertexIndices[currentpos]);

//    ++currentpos;
//    reloadVertices();
//}

void Viewer::handleDeformRightSide(vector<pair<int, vec2> > &shiftmap, vector<bool> &bfixed)
{
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<vec3> new_raw_vertices = raw_vertices[currentpos];

    vector<vec3> targetvt;
    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
    vector<int> count;
    count.resize(new_raw_vertices.size(), 0);

    for(int i=0; i<shiftmap.size(); i++)
    {
        int id = shiftmap[i].first;
        vec2 pt = shiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

        std::array<GLdouble, 3> screen_coords;
        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
        gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
        count[id]++;
    }

    vector<int> handleid;
    vector<vec3> handlevt;

    for(int i=0; i<sidefixid.size(); i++)
    {
        int id = sidefixid[i];
        handleid.push_back(id);
        handlevt.push_back(new_raw_vertices[id]);
        count[id] = -1;
    }

    for(int i=0; i<count.size(); i++)
    {
        if(count[i]<=0) continue;

        {
            handleid.push_back(i);
            vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);
            handlevt.push_back(tv);
        }
    }


    lapdeform->handle_deform_right(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
    //lapdeform->save_deform_info2(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);
    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);

    ++currentpos;
    cout<<currentpos<<" "<<raw_vertices.size()<<endl;
    reloadVertices();
}

void Viewer::handleDeformLeftSide(vector<pair<int, vec2> > &shiftmap, vector<bool> &bfixed)
{
    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;
    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());

    vector<vec3> new_raw_vertices = raw_vertices[currentpos];

    vector<vec3> targetvt;
    targetvt.resize(new_raw_vertices.size(), vec3(0, 0, 0));
    vector<int> count;
    count.resize(new_raw_vertices.size(), 0);

    for(int i=0; i<shiftmap.size(); i++)
    {
        int id = shiftmap[i].first;
        vec2 pt = shiftmap[i].second;

        float x, y, z;
        x = y = z = 0;

        std::array<GLdouble, 3> screen_coords;
        gluProject(new_raw_vertices[id].x(), new_raw_vertices[id].y(), new_raw_vertices[id].z(),
            modelview.data(), projection.data(), viewport.data(),
            screen_coords.data(), screen_coords.data() + 1, screen_coords.data() + 2);

        GLdouble newobj_vertex[3];
        gluUnProject(pt.x(), VIEW_SIZE-pt.y(), screen_coords[2], modelview.data(),
                  projection.data(), viewport.data(), &newobj_vertex[0], &newobj_vertex[1], &newobj_vertex[2]);

        x = newobj_vertex[0], y = newobj_vertex[1], z = newobj_vertex[2];
        targetvt[id][0]+=x;  targetvt[id][1]+=y;  targetvt[id][2]+=z;
        count[id]++;
    }

    vector<int> handleid;
    vector<vec3> handlevt;

    for(int i=0; i<sidefixid.size(); i++)
    {
        int id = sidefixid[i];
        handleid.push_back(id);
        handlevt.push_back(new_raw_vertices[id]);
        count[id] = -1;
    }

    for(int i=0; i<count.size(); i++)
    {
        if(count[i]<=0) continue;

        {
            handleid.push_back(i);
            vec3 tv(0, 0, 0);
            tv[0] = targetvt[i][0]/(1.0f*count[i]);
             tv[1] = targetvt[i][1]/(1.0f*count[i]);
              tv[2] = targetvt[i][2]/(1.0f*count[i]);
            handlevt.push_back(tv);
        }
    }


    lapdeform->handle_deform_left(new_raw_vertices, vertexIndices[currentpos], handleid, handlevt, bfixed);

    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);

    ++currentpos;
    cout<<currentpos<<" "<<raw_vertices.size()<<endl;
    reloadVertices();
}

void Viewer::locateLineArea(int left, int right, int up, int down, vector<QPoint>& line, vector<int>& lpt, vector<float>& lwt)
{
    int np = line.size();

    for (int i = 0; i + 1 < np; ++i){
        QPoint a = line[i];
        QPoint b = line[i+1];
        // Too Narrow
        float draw_dis = (a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y());
        if (draw_dis < 8 && i + 2 != np){
            line.erase(line.begin()+i+1, line.begin()+i+2);
            --i;
            --np;
        }
        // Too Wide
        else if (draw_dis > 40){
            line.insert(line.begin()+i+1, QPoint((a.x() + b.x())/2, (a.y() + b.y())/2));
            --i;
            ++np;
        }
    }



    int nq = raw_vertices[currentpos].size();
    int n_face = vertexIndices[currentpos].size();

    cout<<nq<<" fuck "<<n_face<<" "<<np<<endl;

    int width_extend = (right-left)/4 < 50 ? 50 : (right-left)/4;
    int height_extend = (down-up)/4 < 50 ? 50 : (down-up)/4;
    left < width_extend? left = 0: left -=width_extend;
    up < height_extend? up = 0: up -=height_extend;
    right > VIEW_SIZE - width_extend? right = VIEW_SIZE: right += width_extend;
    down > VIEW_SIZE - height_extend? down = VIEW_SIZE: down += height_extend;

    GLdouble current_x, current_y, current_z;

    std::array<GLdouble, 16> projection;
    std::array<GLdouble, 16> modelview;
    std::array<GLint, 4> viewport;

    glGetDoublev(GL_PROJECTION_MATRIX, projection.data());
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview.data());
    glGetIntegerv(GL_VIEWPORT, viewport.data());
    const vec3 viewdir (camera()->viewDirection().x, camera()->viewDirection().y, camera()->viewDirection().z);
    vec3 campos(camera()->position().x, camera()->position().y, camera()->position().z);

    bdeformflag.clear();
    bdeformflag.resize(nq, false);

    vector<bool> bhandle;
    bhandle.resize(nq, false);
    vector<float> ldis;
    ldis.resize(nq, 0);

    float maxdist, mindist;
    maxdist = -1e+9; mindist = 1e+9;

    float dist_thr = 0;

    vector< vector<float>> dists;
    dists.resize(line.size());

    for(int k=0; k<line.size(); k++)
    {
        dists[k].resize(nq, -1);
        float dist = 1e+5; int id = -1;
        for (int i = 0; i < nq; ++i){
            float cur_angle = QVector3D::dotProduct(viewdir, raw_point_normals[currentpos][i]);
            if(cur_angle < 0){ // the front points
                gluProject(raw_vertices[currentpos][i].x(), raw_vertices[currentpos][i].y(), raw_vertices[currentpos][i].z(),
                    modelview.data(), projection.data(), viewport.data(),
                    &current_x, &current_y, &current_z);
                current_y = VIEW_SIZE - current_y;
                int int_current_x = current_x, int_current_y = current_y;
                if (int_current_x >= left&& int_current_x <= right && int_current_y >= up && int_current_y <= down) {
                     float dis = distPoint(line[k].x(), line[k].y(), int_current_x, VIEW_SIZE-int_current_y);
                     dists[k][i] = dis;
                     if(dis<dist)
                     {
                        dist = dis; id = i;
                     }
                }
            }
        }

        if(id>=0)
        {
            float len = 0; int ncount = 0;
            vec3 p = raw_vertices[currentpos][id];
            GLdouble x0, y0, z0, x1, y1, z1;
            gluProject(p.x(), p.y(), p.z(),modelview.data(), projection.data(), viewport.data(), &x0, &y0, &z0);

            for(auto it = neighbors[currentpos][id].begin(); it!=neighbors[currentpos][id].end(); it++)
            {
               vec3 q = raw_vertices[currentpos][*it];
               gluProject(q.x(), q.y(), q.z(),modelview.data(), projection.data(), viewport.data(), &x1, &y1, &z1);
               len+=sqrtf((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)); ncount++;
            }

            if(ncount>0) len/=ncount;
            if(len>dist_thr) dist_thr = len ;
        }
    }

    for(int i=0; i<nq; i++)
    {
        float dist = 1e+5;
        for(int k=0; k<line.size(); k++)
        {
            if(dists[k][i]>=0&&dists[k][i]<dist)
                dist = dists[k][i];
        }
        if(dist<dist_thr*1.8)
            bdeformflag[i] = true;
    }

    for(int k=0; k<line.size(); k++)
    {
        float dist = 1e+5; int id = -1;
        for(int i=0; i<nq; i++)
        {
            if(!bdeformflag[i]) continue;
            if(dists[k][i]<dist&&dists[k][i]>=0)
            {
                dist = dists[k][i]; id = i;
            }
        }
        if(id>=0) bhandle[id] = true;
    }
/*
    highlight_vertexIndices.clear();
    highlight_vertexIndices.resize(nq, false);
    */
    lpt.clear(); lwt.clear();
    for(int i=0; i<nq; ++i)
    {
        if(bhandle[i])
        {
            /*highlight_vertexIndices[i] = true;*/
            lpt.push_back(i); lwt.push_back(1.0);
        }
        if(bdeformflag[i])
        {
           // highlight_vertexIndices[i] = true;
        }
    }

    update();
}

void Viewer::localdeform()
{
   // smoothSampleArea();

    vector<QVector3D> new_raw_vertices = raw_vertices[currentpos];

    int nvert = new_raw_vertices.size();
    // prepare fflag
    vector<int> fflag; fflag.resize(nvert, -1);

   // bdeformflag = bflag;

    for(int i=0; i<bdeformflag.size(); i++)
    {
        if(bdeformflag[i])
        {
            fflag[i] = 0;
            bool bb = false;
            for(auto it=neighbors[currentpos][i].begin(); it!=neighbors[currentpos][i].end(); it++)
            {
                if(!bdeformflag[*it]) bb = true;
            }
            if(bb) fflag[i] = 1;
        }
    }

    cout<<"start deform"<<endl;
    lapdeform->localDeform(shift_map, vertexIndices[currentpos], new_raw_vertices, neighbors[currentpos], fflag);
    cout<<"finish deform"<<endl;


    for(int iter= 0; iter<1; iter++)
    {
        // smoothing mesh
        for(int i=0; i<new_raw_vertices.size(); i++)
        {
            if(fflag[i]>=0)
            {
                vec3 vt(0, 0, 0); int ncount = 0;
                for(std::set<int>::iterator it = neighbors[currentpos][i].begin();
                    it!=neighbors[currentpos][i].end(); ++it)
                {
                    int jd = *it;
                    vec3 v = new_raw_vertices[jd];
                    for(int k=0; k<3; k++)
                        vt[k]+=v[k];
                    ncount++;
                }
                if(ncount>0)
                {
                    for(int k=0; k<3; k++)
                         new_raw_vertices[i][k]=vt[k]/ncount;
                }
            }
        }
    }

    for(int iter= 0; iter<3; iter++)
    {
        // smoothing mesh
        for(int i=0; i<new_raw_vertices.size(); i++)
        {
            if(fflag[i]>0)
            {
                vec3 vt(0, 0, 0); int ncount = 0;
                for(std::set<int>::iterator it = neighbors[currentpos][i].begin();
                    it!=neighbors[currentpos][i].end(); ++it)
                {
                    int jd = *it;
                    vec3 v = new_raw_vertices[jd];
                    for(int k=0; k<3; k++)
                        vt[k]+=v[k];
                    ncount++;
                }
                if(ncount>0)
                {
                    for(int k=0; k<3; k++)
                         new_raw_vertices[i][k]=vt[k]/ncount;
                }
            }
        }
    }



    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    --currentpos;
    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);


    ++currentpos;
    reloadVertices(); //correct the point normals
}

void Viewer::localROIDeform()
{
    vector<QVector3D> new_raw_vertices = raw_vertices[currentpos];
    int nvert = new_raw_vertices.size();

    vector<bool> bdeformflag;
    bdeformflag.resize(nvert, false);

    for(int i=0; i<nvert; i++)
    {
        if(hightlight_selection_rawindex[i]||hightlight_deformneighbor_rawindex[i])
            bdeformflag[i] = true;
    }
    // prepare fflag
    vector<int> fflag; fflag.resize(nvert, -1);

    for(int i=0; i<bdeformflag.size(); i++)
    {
        if(bdeformflag[i])
        {
            fflag[i] = 0;
            bool bb = false;
            for(auto it=neighbors[currentpos][i].begin(); it!=neighbors[currentpos][i].end(); it++)
            {
                if(!bdeformflag[*it]) bb = true;
            }
            if(bb) fflag[i] = 1;
        }
    }

    for(int i=0; i<nvert; i++)
    {
        bool bb = false;
        if(fflag[i]=-1)
        {
            for(auto it=neighbors[currentpos][i].begin(); it!=neighbors[currentpos][i].end(); it++)
            {
                if(fflag[*it]==1||fflag[*it]==0) bb = true;
            }
            if(bb) fflag[i] = 1;
        }
    }


    lapdeform->localDeform(shift_map, vertexIndices[currentpos], new_raw_vertices, neighbors[currentpos], fflag);

    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);

    ++currentpos;
    reloadVertices(); //correct the point normals
}

void Viewer::deform(bool partial, bool shift, int level, const std::vector<int> & selectedNeighbors){

    vector<QVector3D> new_raw_vertices = raw_vertices[currentpos];

    lapdeform->lapExaggeration(shift_map, vertexIndices[currentpos], new_raw_vertices, neighbors[currentpos], LAP_NO_DEBUG, partial, shift, level, selectedNeighbors);


    vector<int> new_vertexIndices = vertexIndices[currentpos];
    vector<set<int>> new_neighbors = neighbors[currentpos];

    --currentpos;
    removeAfterCurrentpos();

    raw_vertices.push_back(new_raw_vertices);
    vertexIndices.push_back(new_vertexIndices);
    neighbors.push_back(new_neighbors);


    ++currentpos;
    reloadVertices(); //correct the point normals
}

void Viewer::removeAfterCurrentpos(){
    if (currentpos+1 == vertexIndices.size()) return;
    vertexIndices.erase(vertexIndices.begin()+currentpos+1, vertexIndices.end());
    raw_vertices.erase(raw_vertices.begin()+currentpos+1, raw_vertices.end());
    vertices.erase(vertices.begin()+currentpos+1, vertices.end());
    raw_point_normals.erase(raw_point_normals.begin()+currentpos+1, raw_point_normals.end());
    point_normals.erase(point_normals.begin()+currentpos+1, point_normals.end());
    poly_normals.erase(poly_normals.begin()+currentpos+1, poly_normals.end());
    neighbors.erase(neighbors.begin()+currentpos+1, neighbors.end());
}

void Viewer::draw() {

    int nvert = vertices[currentpos].size();
    int nq = raw_vertices[currentpos].size();

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    // wire
    if (wire){
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

        glColor3f(0.40f, 0.40f, 0.40f);
        glBegin(GL_LINES);
        for (int i=0; i<nvert; i+=3){
            glVertex3f(vertices[currentpos][i].x(), vertices[currentpos][i].y(), vertices[currentpos][i].z());
            glVertex3f(vertices[currentpos][i+1].x(), vertices[currentpos][i+1].y(), vertices[currentpos][i+1].z());
            glVertex3f(vertices[currentpos][i+1].x(), vertices[currentpos][i+1].y(), vertices[currentpos][i+1].z());
            glVertex3f(vertices[currentpos][i+2].x(), vertices[currentpos][i+2].y(), vertices[currentpos][i+2].z());
            glVertex3f(vertices[currentpos][i+2].x(), vertices[currentpos][i+2].y(), vertices[currentpos][i+2].z());
            glVertex3f(vertices[currentpos][i].x(), vertices[currentpos][i].y(), vertices[currentpos][i].z());
        }
        glEnd();
        glDisable(GL_BLEND);
        glDisable(GL_LINE_SMOOTH);
    }

    float hr, hg, hb;
    float r, g, b;

    hr = 0.8f;  hg = 0.4f;  hb = 0.4f;
    r = 0.7f; g = 0.7f; b = 0.7f;

    if (smooth){
        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

        glBegin(GL_TRIANGLES);
        glColor3f(r, g, b);
        for (int i=0; i<nvert; i+=3){
            glNormal3f(point_normals[currentpos][i].x(), point_normals[currentpos][i].y(), point_normals[currentpos][i].z());
            glVertex3f(vertices[currentpos][i].x(), vertices[currentpos][i].y(), vertices[currentpos][i].z());
            glNormal3f(point_normals[currentpos][i+1].x(), point_normals[currentpos][i+1].y(), point_normals[currentpos][i+1].z());
            glVertex3f(vertices[currentpos][i+1].x(), vertices[currentpos][i+1].y(), vertices[currentpos][i+1].z());
            glNormal3f(point_normals[currentpos][i+2].x(), point_normals[currentpos][i+2].y(), point_normals[currentpos][i+2].z());
            glVertex3f(vertices[currentpos][i+2].x(), vertices[currentpos][i+2].y(), vertices[currentpos][i+2].z());
        }
        glEnd();

        if(isROISelection)
        {
            glPointSize(3);
            glColor3f(0.6, 0.6, 0.6);
            glBegin(GL_POINTS);
            for (int i=0; i<nq; ++i){
                if (!hightlight_deformneighbor_rawindex[i]){
                    glVertex3f(raw_vertices[currentpos][i].x(),raw_vertices[currentpos][i].y(),raw_vertices[currentpos][i].z());
                }
            }
            glEnd();
        }

        if(isHandleSelection)
        {
            glPointSize(8);
            glColor3f(hr, hg, hb);
            glBegin(GL_POINTS);
            for (int i=0; i<nq; ++i){
                if (hightlight_selection_rawindex[i]){
                    glVertex3f(raw_vertices[currentpos][i].x(),raw_vertices[currentpos][i].y(),raw_vertices[currentpos][i].z());
                }
            }
            glEnd();
        }

        glDisable(GL_POLYGON_SMOOTH);
        glColor3f(r, g, b);
    }


    if (texture){

        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

        // Enable GL textures
        glEnable(GL_TEXTURE_2D);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

        glBegin(GL_TRIANGLES);
        for (int i=0; i<nvert; ++i){
            glTexCoord2f(vts[vts_label[vertexIndices[currentpos][i]-1]].first, vts[vts_label[vertexIndices[currentpos][i]-1]].second);
            glNormal3f(point_normals[currentpos][i].x(), point_normals[currentpos][i].y(), point_normals[currentpos][i].z());
            glVertex3f(vertices[currentpos][i].x(), vertices[currentpos][i].y(), vertices[currentpos][i].z());
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);

        if(isROISelection)
        {
            glPointSize(3);
            glColor3f(0.6, 0.6, 0.6);
            glBegin(GL_POINTS);
            for (int i=0; i<nq; ++i){
                if (!hightlight_deformneighbor_rawindex[i]){
                    glVertex3f(raw_vertices[currentpos][i].x(),raw_vertices[currentpos][i].y(),raw_vertices[currentpos][i].z());
                }
            }
            glEnd();
        }

        if(isHandleSelection)
        {
            glPointSize(8);
            glColor3f(hr, hg, hb);
            glBegin(GL_POINTS);
            for (int i=0; i<nq; ++i){
                if (hightlight_selection_rawindex[i]){
                    glVertex3f(raw_vertices[currentpos][i].x(),raw_vertices[currentpos][i].y(),raw_vertices[currentpos][i].z());
                }
            }
            glEnd();
        }

        glDisable(GL_POLYGON_SMOOTH);
        glColor3f(r, g, b);
    }

    updatePlayers();
    glDisable(GL_LIGHTING);
    displayPlayers();
    glEnable(GL_LIGHTING);
}

void Viewer::setEnvironment(unsigned int bgr, unsigned int bgg, unsigned int bgb){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    QColor q(bgr,bgg,bgb);
    setBackgroundColor(q);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);


   // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    fLightAmbient[0]  = fLightAmbient[1]  = fLightAmbient[2]  = 0.5f;
    fLightDiffuse[0]  = fLightDiffuse[1]  = fLightDiffuse[2]  =0.5f;
    fLightSpecular[0] = fLightSpecular[1] = fLightSpecular[2] = 0.5f;
    fLightAmbient[3]  = fLightDiffuse[3]  = fLightSpecular[3] = 1.0f;
    fLightShininess = 0.5f;
    fLightPosition[0] = 10.0;	fLightPosition[1] = 10.0;
    fLightPosition[2] = 20.0;	fLightPosition[3] = 0.5;

    SetLightParamsV(GL_LIGHT0,fLightAmbient,fLightDiffuse,fLightSpecular,fLightShininess);
    glLightfv(GL_LIGHT0, GL_POSITION, fLightPosition);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
    glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0);

    setMatPrameter("Pewter");
    SetMaterialParamsV(GL_FRONT,fMatAmbientFront,fMatDiffuseFront,fMatSpecularFront);
    for (int i=0; i<4; ++i)		fMatAmbientBack[i] = 1.0f * fMatAmbientFront[i];
    for (int i=0; i<4; ++i)		fMatDiffuseBack[i] = 1.0f * fMatDiffuseFront[i];
    for (int i=0; i<4; ++i)		fMatSpecularBack[i] = 1.0f * fMatSpecularFront[i];
    SetMaterialParamsV(GL_BACK,fMatAmbientBack,fMatDiffuseBack,fMatSpecularBack);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, fMatShininess);

    camera()->setSceneCenter(sceneCenter());
    camera()->loadModelViewMatrix();

    const GLdouble modelview[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, -0.2, -4.7, 1};
    camera()->setFromModelViewMatrix(modelview);

}

void Viewer::MoveToFront(){

    camera()->setSceneCenter(sceneCenter());
    camera()->loadModelViewMatrix();

    const GLdouble modelview[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, -0.2, -4.7, 1};
    camera()->setFromModelViewMatrix(modelview);
    update();
}

void Viewer::MoveToSide(){

    camera()->setSceneCenter(sceneCenter());
    camera()->loadModelViewMatrix();

    const GLdouble modelview[16] = {0.03187, -0.01831, 0.999, 0,
                                    -0.068579, 0.997, 0.02046, 0,
                                    -0.997, -0.06918, 0.03054, 0,
                                    -0.07098, -0.2, -4.7, 1
                                   };
    camera()->setFromModelViewMatrix(modelview);
    update();

}

void Viewer::MoveToLeftSide(){

    camera()->setSceneCenter(sceneCenter());
    camera()->loadModelViewMatrix();

    const GLdouble modelview[16] = {0, 0, -1, 0,
                                    0, 1, 0, 0,
                                    1, 0, 0, 0,
                                    0, -0.2, -4.7, 1
                                   };
    camera()->setFromModelViewMatrix(modelview);
    update();

}

void Viewer::MoveToRightSide(){

    camera()->setSceneCenter(sceneCenter());
    camera()->loadModelViewMatrix();

    const GLdouble modelview[16] = {0, 0, 1, 0,
                                    0, 1, 0, 0,
                                    -1, 0, 0, 0,
                                    0, -0.223663, -4.76, 1
                                   };
    camera()->setFromModelViewMatrix(modelview);
    update();

}

void CameraPathPlayer::checkIfGrabsMouse(int x, int y, const Camera* const){
  setGrabsMouse((x < 80) && (y<yPos()) && ((yPos()-y) < 16));
}

void Viewer::displayPlayers(){
    for (int i=0; i<nbPlayers_; ++i){
        CameraPathPlayer* cpp = player_[i];
        if (cpp){
            QString s;
            if (cpp->grabsMouse()){
                glColor3f(1,1,1);
                if (camera()->keyFrameInterpolator(i)->numberOfKeyFrames() > 1)
                    s = "Play path F" + QString::number(i);
                else
                    s = "Restore pos F" + QString::number(i);
            }
            else {
                glColor3f(0.6f, 0.6f, 0.6f);
                if (camera()->keyFrameInterpolator(i)->numberOfKeyFrames() > 1)
                    s = "Path F" + QString::number(i);
                else
                    s = "Pos F" + QString::number(i);
            }
             drawText(10, cpp->yPos()-3, s);
        }
    }
}

void Viewer::updatePlayers() {
    for (int i=0; i<nbPlayers_; ++i){
        // Check if CameraPathPlayer is still valid
        if ((player_[i]) && (!camera()->keyFrameInterpolator(i))){
            delete player_[i];
            player_[i] = NULL;
        }
        // Or add it if needed
        if ((camera()->keyFrameInterpolator(i)) && (!player_[i])){
            player_[i] = new CameraPathPlayer(i);
        }
    }
}

void Viewer::setMatPrameter(const char* name) {
  GLfloat ambient[] = {0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat diffuse[] = {0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat specular[] = {0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat shininess[] = {0.0f};

  if (name == "Silver") {
    Set4Array(ambient, 0.19225f, 0.19225f, 0.19225f, 1.0f);
    Set4Array(diffuse, 0.50754f, 0.50754f, 0.50754f, 1.0f);
    Set4Array(specular, 0.208273f, 0.208273f, 0.208273f, 1.0f);
    shininess[0] = 51.2f;
  } else if (name == "Gold") {
    Set4Array(ambient, 0.24725f, 0.1995f, 0.0745f, 1.0f);
    Set4Array(diffuse, 0.75164f, 0.60648f, 0.22648f, 1.0f);
    Set4Array(specular, 0.628281f, 0.555802f, 0.366065f, 1.0f);
    shininess[0] = 51.2f;
  } else if (name == "Jade") {
    Set4Array(ambient, 0.135f, 0.2225f, 0.1575f, 0.95f);
    Set4Array(diffuse, 0.54f, 0.89f, 0.63f, 0.95f);
    Set4Array(specular, 0.316228f, 0.316228f, 0.316228f, 0.95f);
    shininess[0] = 12.8f;
  } else if (name == "Light blue") {
    Set4Array(ambient, 0.0f, 0.5f, 0.75f, 1.0f);
    Set4Array(diffuse, 0.0f, 0.5f, 1.0f, 1.0f);
    Set4Array(specular, 0.75f, 0.75f, 0.75f, 1.0f);
    shininess[0] = 64.0f;
  } else if (name == "Cyan blue") {
    Set4Array(ambient, 0.207843f, 0.501961f, 0.486275f, 1.0f);
    Set4Array(diffuse, 0.513726f, 0.921569f, 0.85098f, 1.0f);
    Set4Array(specular, 0.309804f, 0.309804f, 0.309804f, 1.0f);
    shininess[0] = 64.0f;
  } else if (name == "Skin") {
    // Ambient
    ambient[0] = 0.83f;
    ambient[1] = 0.69f;
    ambient[2] = 0.588f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.83f;
    diffuse[1] = 0.69f;
    diffuse[2] = 0.588f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.25f;
    specular[1] = 0.25f;
    specular[2] = 0.25f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 64.0f;
  } else if (name == "Emerald") {
    // Ambient
    ambient[0] = 0.0215f;
    ambient[1] = 0.1745f;
    ambient[2] = 0.0215f;
    ambient[3] = 0.55f;
    // Diffuse
    diffuse[0] = 0.07568f;
    diffuse[1] = 0.61424f;
    diffuse[2] = 0.07568f;
    diffuse[3] = 0.55f;
    // Specular
    specular[0] = 0.633f;
    specular[1] = 0.727811f;
    specular[2] = 0.633f;
    specular[3] = 0.55f;
    // Shininess
    shininess[0] = 76.8f;
  } else if (name == "Polished silver") {
    // Ambient
    ambient[0] = 0.23125f;
    ambient[1] = 0.23125f;
    ambient[2] = 0.23125f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.2775f;
    diffuse[1] = 0.2775f;
    diffuse[2] = 0.2775f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.773911f;
    specular[1] = 0.773911f;
    specular[2] = 0.773911f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 89.6f;
  } else if (name == "Chrome") {
    // Ambient
    ambient[0] = 0.25f;
    ambient[1] = 0.25f;
    ambient[2] = 0.25f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.4f;
    diffuse[1] = 0.4f;
    diffuse[2] = 0.4f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.774597f;
    specular[1] = 0.774597f;
    specular[2] = 0.774597f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 76.8f;
  }

  else if (name == "Copper") {
    // Ambient
    ambient[0] = 0.19125f;
    ambient[1] = 0.0735f;
    ambient[2] = 0.0225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.7038f;
    diffuse[1] = 0.27048f;
    diffuse[2] = 0.0828f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.256777f;
    specular[1] = 0.137622f;
    specular[2] = 0.086014f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 12.8f;
  }

  else if (name == "Polished gold") {
    // Ambient
    ambient[0] = 0.24725f;
    ambient[1] = 0.2245f;
    ambient[2] = 0.0645f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.34615f;
    diffuse[1] = 0.3143f;
    diffuse[2] = 0.0903f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.797357f;
    specular[1] = 0.723991f;
    specular[2] = 0.208006f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 83.2f;
  }

  else if (name == "Pewter") {
    // Ambient
    ambient[0] = 0.105882f;
    ambient[1] = 0.058824f;
    ambient[2] = 0.113725f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.427451f;
    diffuse[1] = 0.470588f;
    diffuse[2] = 0.541176f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.333333f;
    specular[1] = 0.333333f;
    specular[2] = 0.521569f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 9.84615f;
  }

  else if (name == "Obsidian") {
    // Ambient
    ambient[0] = 0.05375f;
    ambient[1] = 0.05f;
    ambient[2] = 0.06625f;
    ambient[3] = 0.82f;
    // Diffuse
    diffuse[0] = 0.18275f;
    diffuse[1] = 0.17f;
    diffuse[2] = 0.22525f;
    diffuse[3] = 0.82f;
    // Specular
    specular[0] = 0.332741f;
    specular[1] = 0.328634f;
    specular[2] = 0.346435f;
    specular[3] = 0.82f;
    // Shininess
    shininess[0] = 38.4f;
  }

  else if (name == "Black plastic") {
    // Ambient
    ambient[0] = 0.0f;
    ambient[1] = 0.0f;
    ambient[2] = 0.0f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.01f;
    diffuse[1] = 0.01f;
    diffuse[2] = 0.01f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.5f;
    specular[1] = 0.5f;
    specular[2] = 0.5f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 32.0f;
  }

  else if (name == "Polished bronze") {
    // Ambient
    ambient[0] = 0.25f;
    ambient[1] = 0.148f;
    ambient[2] = 0.006475f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.4f;
    diffuse[1] = 0.2368f;
    diffuse[2] = 0.1036f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.774597f;
    specular[1] = 0.458561f;
    specular[2] = 0.200621f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 76.8f;
  }

  else if (name == "Polished copper") {
    // Ambient
    ambient[0] = 0.2295f;
    ambient[1] = 0.08825f;
    ambient[2] = 0.0275f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.5508f;
    diffuse[1] = 0.2118f;
    diffuse[2] = 0.066f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.580594f;
    specular[1] = 0.223257f;
    specular[2] = 0.0695701f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 51.2f;
  }

  else if (name == "Pearl") {
    // Ambient
    ambient[0] = 0.25f;
    ambient[1] = 0.20725f;
    ambient[2] = 0.20725f;
    ambient[3] = 0.922f;
    // Diffuse
    diffuse[0] = 1.0f;
    diffuse[1] = 0.829f;
    diffuse[2] = 0.829f;
    diffuse[3] = 0.922f;
    // Specular
    specular[0] = 0.296648f;
    specular[1] = 0.296648f;
    specular[2] = 0.296648f;
    specular[3] = 0.922f;
    // Shininess
    shininess[0] = 11.264f;
  }

  else if (name == "Ruby") {
    // Ambient
    ambient[0] = 0.1745f;
    ambient[1] = 0.01175f;
    ambient[2] = 0.01175f;
    ambient[3] = 0.55f;
    // Diffuse
    diffuse[0] = 0.61424f;
    diffuse[1] = 0.04136f;
    diffuse[2] = 0.04136f;
    diffuse[3] = 0.55f;
    // Specular
    specular[0] = 0.727811f;
    specular[1] = 0.626959f;
    specular[2] = 0.626959f;
    specular[3] = 0.55f;
    // Shininess
    shininess[0] = 76.8f;
  }

  else if (name == "Turquoise") {
    // Ambient
    ambient[0] = 0.1f;
    ambient[1] = 0.18725f;
    ambient[2] = 0.1745f;
    ambient[3] = 0.8f;
    // Diffuse
    diffuse[0] = 0.396f;
    diffuse[1] = 0.74151f;
    diffuse[2] = 0.69102f;
    diffuse[3] = 0.8f;
    // Specular
    specular[0] = 0.297254f;
    specular[1] = 0.30829f;
    specular[2] = 0.306678f;
    specular[3] = 0.8f;
    // Shininess
    shininess[0] = 12.8f;
  }

  else if (name == "Brass") {
    // Ambient
    ambient[0] = 0.329412f;
    ambient[1] = 0.223529f;
    ambient[2] = 0.027451f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.780392f;
    diffuse[1] = 0.268627f;
    diffuse[2] = 0.113725f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.992157f;
    specular[1] = 0.741176f;
    specular[2] = 0.807843f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 27.8974f;
  }

  else if (name == "Brass") {
    // Ambient
    ambient[0] = 0.329412f;
    ambient[1] = 0.223529f;
    ambient[2] = 0.027451f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.780392f;
    diffuse[1] = 0.268627f;
    diffuse[2] = 0.113725f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.992157f;
    specular[1] = 0.741176f;
    specular[2] = 0.807843f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 27.8974f;
  }

  for (int i = 0; i < 4; i++) {
    fMatAmbientFront[i] = ambient[i];
    fMatDiffuseFront[i] = diffuse[i];
    fMatSpecularFront[i] = specular[i];
  }
  fMatShininess = shininess[0];
}
