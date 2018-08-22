#ifndef DEFORM_H
#define DEFORM_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>

#include <vector>
#include <QVector3D>
#include <map>
#include <set>
#include <unordered_set>

#include <igl/grad.h>
#include <igl/cotmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/writeOBJ.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/writeDMAT.h>
#include <igl/doublearea.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//#include <cuda_runtime.h>
//#include <cusparse_v2.h>
//#include <cusolverSp.h>
#include <igl/harmonic.h>
#include <fstream>

//#include "unitilies.h"
#include "folder.h"

using namespace std;

typedef QVector3D vec3;

#define LAP_DEBUG true
#define LAP_NO_DEBUG false
#define LAPMAT_DIM 80134
#define NVERT 11510
#define NFIXED 669
#define NFIXEDLEFT 1645
#define NFIXEDRIGHT 1644

class LapDeform
{
public:
    LapDeform(){
         initblapfixed(std::string(localfolder()+"\\models\\lapDeform\\idfixedlap3.in").c_str());
         initHandleLeft(std::string(localfolder()+"\\models\\lapDeform\\handleidleft.in").c_str());
         initLapmatLeft(std::string(localfolder()+"\\models\\lapDeform\\lapmat.tri").c_str());
         initHandleRight(std::string(localfolder()+"\\models\\lapDeform\\handleidright.in").c_str());
         initLapmatRight(std::string(localfolder()+"\\models\\lapDeform\\lapmat.tri").c_str());
         initHandle(std::string(localfolder()+"\\models\\lapDeform\\handleid.in").c_str());
         initLapmat(std::string(localfolder()+"\\models\\lapDeform\\lapmat.tri").c_str());
         printf("Lap Done\n");
    }
    ~LapDeform(){}

    void initLapmat(const char* file);
    void initblapfixed(const char* file);
    void initHandle(const char* file);
    void initHandleLeft(const char* file);
    void initLapmatLeft(const char* file);
    void initHandleRight(const char* file);
    void initLapmatRight(const char* file);

    void add_neighbors(const std::vector<QVector3D> &raw_vertices,
                       std::vector<int> & vertex_candidates,
                       std::vector<int> & kp,
                       std::vector<QVector3D> & kvt,
                       const std::vector<std::set<int> > & neighbors,
                       int currrent_index,
                       int iteration);

    void lapExaggeration(const std::map<int, QVector3D> & shift_map,
                         const std::vector<int> & vertexIndices,
                         std::vector<QVector3D> & raw_vertices,
                         const std::vector<std::set<int> > & neighbors,
                         bool debug_flag,
                         bool partial_flag,
                         bool shift_flag,
                         int level,
                         const std::vector<int> & selectedNeighbors);

    void localDeform(map<int, vec3>& shift_map, vector<int>& vertexIndices,
                         vector<vec3>& raw_vertices, vector<std::set<int>>& neighbors, vector<int>& fflag);

    void handle_deform(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                           std::vector<vec3>& handlePositions, vector<bool>& bfixed);

    void handle_deform_right(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                             std::vector<vec3>& handlePositions, vector<bool>& bfixed);

    void handle_deform_left(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                             std::vector<vec3>& handlePositions, vector<bool>& bfixed);


    void handle_deform_all(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                           std::vector<vec3>& handlePositions, vector<bool>& bfixed);

    void handle_deform_cuda(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                           std::vector<vec3>& handlePositions);

    void save_deform_info(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                          std::vector<vec3>& handlePositions, vector<bool>& bfixed);

//    void save_deform_info2(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
//                          std::vector<vec3>& handlePositions, vector<bool>& bfixed);

    void lapDeform_cuda(const std::map<int, QVector3D> & shift_map,
                         const std::vector<int> & vertexIndices,
                         std::vector<QVector3D> & raw_vertices,
                         const std::vector<std::set<int> > & neighbors,
                         bool debug_flag,
                         bool partial_flag,
                         bool shift_flag);
    vector<int> vertex_candidates;

private:
    Eigen::SparseMatrix<double> lapmat, lapmatleft, coffemat, coffematleft, lapmatright, coffematright;
    vector<float> handle_weight, handle_weight_left, handle_weight_right;
    vector<int> handle_id, handle_id_left, handle_id_right;
    vector<int> lap_i, lap_j;
    vector<double> lap_v;
    vector<bool> blapfixed;
    Eigen::SparseMatrix<double> Amat, AmatLeft, AmatRight;
    Eigen::SparseMatrix<double> mat_a, mat_a_left, mat_a_right;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> lapsolver, lapsolver_left, lapsolver_right;

    unordered_set<int> visited_neighbors;
    int neighbor_range;

};


//partial_flag: all points may move = false, only move neighbors = true
//shift_flag: shiftmode(silhouette) = true, dragmode = false;

#endif // DEFORM_H
