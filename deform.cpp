#include "deform.h"

void LapDeform::initLapmat(const char *file)
{
    lap_i.resize(LAPMAT_DIM);
    lap_j.resize(LAPMAT_DIM);
    lap_v.resize(LAPMAT_DIM);
    ifstream lap(file);
    for(int k=0; k<LAPMAT_DIM; k++)
    {
        int i, j;
        double value;
        lap>>i>>j>>value;
        lap_i[k] = i;
        lap_j[k] = j;
        lap_v[k] = value;
    }
    lap.close();

    // init lap matrix
    lapmat.resize(NVERT, NVERT);
    lapmat.reserve(Eigen::VectorXi::Constant(NVERT,10));

    for(int k=0; k<lap_i.size(); k++)
    {
       lapmat.insert(lap_i[k], lap_j[k]) = lap_v[k];
    }

    // init coffemat
    coffemat.resize(NVERT+NFIXED, NVERT);
    coffemat.reserve(Eigen::VectorXi::Constant(NVERT+NFIXED,10));
    for(int k=0; k<lap_i.size(); k++)
    {
       float weight = 1.0;
       int id = lap_i[k];
       if(blapfixed[id]) weight = 3.0;
       coffemat.insert(lap_i[k], lap_j[k]) = weight*lap_v[k];
    }

    for(int k=0; k<handle_id.size(); k++)
    {
        int id = handle_id[k];
        float wt = handle_weight[k];
        coffemat.insert(NVERT+k, id) = wt;
    }

    Amat.resize(coffemat.rows(), coffemat.cols());
    Amat = coffemat;
    mat_a = Amat.transpose() * Amat;
    lapsolver.compute(mat_a);
}

void LapDeform::initLapmatLeft(const char *file)
{
    lap_i.resize(LAPMAT_DIM);
    lap_j.resize(LAPMAT_DIM);
    lap_v.resize(LAPMAT_DIM);
    ifstream lap(file);
    for(int k=0; k<LAPMAT_DIM; k++)
    {
        int i, j;
        double value;
        lap>>i>>j>>value;
        lap_i[k] = i;
        lap_j[k] = j;
        lap_v[k] = value;
    }
    lap.close();

    // init lap matrix
    lapmatleft.resize(NVERT, NVERT);
    lapmatleft.reserve(Eigen::VectorXi::Constant(NVERT,10));

    for(int k=0; k<lap_i.size(); k++)
    {
       lapmatleft.insert(lap_i[k], lap_j[k]) = lap_v[k];
    }

    // init coffemat
    coffematleft.resize(NVERT+NFIXEDLEFT, NVERT);
    coffematleft.reserve(Eigen::VectorXi::Constant(NVERT+NFIXEDLEFT,10));
    for(int k=0; k<lap_i.size(); k++)
    {
       float weight = 1.0;
       int id = lap_i[k];
       if(blapfixed[id]) weight = 3.0;
       coffematleft.insert(lap_i[k], lap_j[k]) = weight*lap_v[k];
    }

    for(int k=0; k<handle_id_left.size(); k++)
    {
        int id = handle_id_left[k];
        float wt = handle_weight_left[k];
        coffematleft.insert(NVERT+k, id) = wt;
    }

    AmatLeft.resize(coffematleft.rows(), coffematleft.cols());
    AmatLeft = coffematleft;
    mat_a_left = AmatLeft.transpose() * AmatLeft;
    lapsolver_left.compute(mat_a_left);
}

void LapDeform::initLapmatRight(const char *file)
{
    lap_i.resize(LAPMAT_DIM);
    lap_j.resize(LAPMAT_DIM);
    lap_v.resize(LAPMAT_DIM);
    ifstream lap(file);
    for(int k=0; k<LAPMAT_DIM; k++)
    {
        int i, j;
        double value;
        lap>>i>>j>>value;
        lap_i[k] = i;
        lap_j[k] = j;
        lap_v[k] = value;
    }
    lap.close();

    // init lap matrix
    lapmatright.resize(NVERT, NVERT);
    lapmatright.reserve(Eigen::VectorXi::Constant(NVERT,10));

    for(int k=0; k<lap_i.size(); k++)
    {
       lapmatright.insert(lap_i[k], lap_j[k]) = lap_v[k];
    }

    // init coffemat
    coffematright.resize(NVERT+NFIXEDLEFT, NVERT);
    coffematright.reserve(Eigen::VectorXi::Constant(NVERT+NFIXEDRIGHT,10));
    for(int k=0; k<lap_i.size(); k++)
    {
       float weight = 1.0;
       int id = lap_i[k];
       if(blapfixed[id]) weight = 3.0;
       coffematright.insert(lap_i[k], lap_j[k]) = weight*lap_v[k];
    }

    for(int k=0; k<handle_id_right.size(); k++)
    {
        int id = handle_id_right[k];
        float wt = handle_weight_right[k];
        coffematright.insert(NVERT+k, id) = wt;
    }

    AmatRight.resize(coffematright.rows(), coffematright.cols());
    AmatRight = coffematright;
    mat_a_right = AmatRight.transpose() * AmatRight;
    lapsolver_right.compute(mat_a_right);
}

void LapDeform::initHandle(const char* file)
{
    handle_weight.clear();
    handle_id.clear();
    handle_weight.resize(NFIXED, 1.0f);
    handle_id.resize(NFIXED, 0);

    ifstream handle(file);
    for(int i=0; i<NFIXED; i++)
    {
        int id, wt;
        handle>>id>>wt;
        handle_id[i] = id;
        handle_weight[i] = wt;
    }
    handle.close();
}

void LapDeform::initHandleLeft(const char* file)
{
    handle_weight_left.clear();
    handle_id_left.clear();
    handle_weight_left.resize(NFIXEDLEFT, 1.0f);
    handle_id_left.resize(NFIXEDLEFT, 0);

    ifstream handle(file);
    for(int i=0; i<NFIXEDLEFT; i++)
    {
        int id, wt;
        handle>>id>>wt;
        handle_id_left[i] = id;
        handle_weight_left[i] = wt;
    }
    handle.close();
}

void LapDeform::initHandleRight(const char* file)
{
    handle_weight_right.clear();
    handle_id_right.clear();
    handle_weight_right.resize(NFIXEDRIGHT, 1.0f);
    handle_id_right.resize(NFIXEDRIGHT, 0);

    ifstream handle(file);
    for(int i=0; i<NFIXEDRIGHT; i++)
    {
        int id, wt;
        handle>>id>>wt;
        handle_id_right[i] = id;
        handle_weight_right[i] = wt;
    }
    handle.close();
}

void LapDeform::initblapfixed(const char *file)
{
    blapfixed.resize(NVERT, false);
    ifstream lapfix(file);
    for(int i=0; i<2682; i++)
    {
        int id; lapfix>>id;
        blapfixed[id] = true;
    }
    lapfix.close();
}

void LapDeform::add_neighbors(const vector<QVector3D> &raw_vertices,
                   vector<int> & vertex_candidates,
                   vector<int> & kp,
                   vector<QVector3D> & kvt,
                   const std::vector<std::set<int> > & neighbors,
                   int currrent_index,
                   int iteration){
//    printf("%d ", iteration);
    if (!iteration) return;
    if (find(vertex_candidates.begin(), vertex_candidates.end(), currrent_index) == vertex_candidates.end()){
        vertex_candidates.push_back(currrent_index);
        // add the last level of vertices to the control point group. for smooth edge
        if (iteration == 1){
            kp.push_back(vertex_candidates.size()-1);
            kvt.push_back(raw_vertices[currrent_index]);
        }
    }

    if (visited_neighbors.find(currrent_index) == visited_neighbors.end()){
        visited_neighbors.insert(currrent_index);
        for (auto it = neighbors[currrent_index].begin(); it != neighbors[currrent_index].end(); ++it){
            add_neighbors(raw_vertices, vertex_candidates, kp, kvt, neighbors, *it, iteration-1);
        }
    } else if (iteration == neighbor_range){
        for (auto it = neighbors[currrent_index].begin(); it != neighbors[currrent_index].end(); ++it){
            add_neighbors(raw_vertices, vertex_candidates, kp, kvt, neighbors, *it, iteration-1);
        }
    }

}

void LapDeform::lapExaggeration(const std::map<int, QVector3D> & shift_map,
                     const std::vector<int> & vertexIndices,
                     std::vector<QVector3D> & raw_vertices,
                     const std::vector<std::set<int> > & neighbors,
                     bool debug_flag,
                     bool partial_flag,
                     bool shift_flag,
                     int level,
                     const std::vector<int> & selectedNeighbors){

    if (debug_flag){
        for(auto it = shift_map.begin(); it != shift_map.end(); ++it) {
          raw_vertices[it->first] = it->second;
        }
        return;
    }

    neighbor_range = shift_flag ? level : 2+level;

//    printf("Start lapExaggeration!\n");

    vector<int> kp; vector<QVector3D> kvt;
    for(auto it = shift_map.begin(); it != shift_map.end(); ++it) {
      kp.push_back(it->first);
      kvt.push_back(it->second);
    }

    int nvert;
    int nface;
    Eigen::MatrixXd v;
    Eigen::MatrixXd f;

    if (partial_flag){
        int nkp = kp.size();
        for (int i = 0; i < nkp; ++i){
            vertex_candidates.push_back(kp[i]);
        }


        //const clock_t start_time = clock(); printf("####Timer Start\n\n");
        //visited_neighbors.clear();
        for (int i = 0; i < nkp; ++i){
            //add_neighbors(raw_vertices, vertex_candidates, kp, kvt, neighbors, kp[i], neighbor_range);
            kp[i] = i; // change to 1 2 3.... order
        }
        //visited_neighbors.clear();
        //printf("Neighbor Total Time: %f\n\nReload:\n", float( clock () - start_time ) /  CLOCKS_PER_SEC);

        int ineighbor = selectedNeighbors.size();
        for (int i = 0; i < ineighbor; ++i) vertex_candidates.push_back(selectedNeighbors[i]);


        nvert = vertex_candidates.size();
        v.resize(nvert, 3);
        int vi = 0;
        for (auto it = vertex_candidates.begin(); it != vertex_candidates.end(); ++it){
            v(vi, 0) = raw_vertices[*it].x();
            v(vi, 1) = raw_vertices[*it].y();
            v(vi, 2) = raw_vertices[*it].z();
            ++vi;
        }
        // select candidate faces

        int temp_nface = vertexIndices.size()/3;
        vector<vector<int> > temp_face;
        for (int i=0; i<temp_nface; i++) {
            int ida = vertexIndices[3*i] - 1;
            int idb = vertexIndices[3*i+1] - 1;
            int idc = vertexIndices[3*i+2] - 1;
            auto ita = std::find(vertex_candidates.begin(), vertex_candidates.end(), ida); if (ita == vertex_candidates.end()) continue;
            auto itb = std::find(vertex_candidates.begin(), vertex_candidates.end(), idb); if (itb == vertex_candidates.end()) continue;
            auto itc = std::find(vertex_candidates.begin(), vertex_candidates.end(), idc); if (itc == vertex_candidates.end()) continue;
            temp_face.push_back(vector<int>{
                                    int(ita - vertex_candidates.begin()),
                                    int(itb - vertex_candidates.begin()),
                                    int(itc - vertex_candidates.begin())});
        }
        nface = temp_face.size();
        f.resize(nface, 3);
        for (int i=0; i<nface; i++) {
            for (int j=0; j<3; j++) {
                f(i, j) = temp_face[i][j];
            }
        }
//        printf("Total %d vertices, %d faces.\n", nvert, nface);
    } else {
        nvert = raw_vertices.size();
        nface = vertexIndices.size()/3;

        v.resize(nvert, 3);

        for (int i=0; i<nvert; ++i) {
            v(i, 0) = raw_vertices[i].x();
            v(i, 1) = raw_vertices[i].y();
            v(i, 2) = raw_vertices[i].z();
        }


        f.resize(nface, 3);

        for (int i=0; i<nface; i++) {
            for (int j=0; j<3; j++) {
                f(i, j) = vertexIndices[3*i+j] - 1;
            }
        }
    }
//    printf("Finish copying data!\n");

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::SparseMatrix<double> G, laplacian_mat;
    igl::grad(v,f,G);
    Eigen::VectorXd dblA;
    igl::doublearea(v,f,dblA);
    const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
    laplacian_mat = G.transpose() * T * G;

    Eigen::MatrixXd lap = laplacian_mat*v;


//    printf("Finish buding lap!\n");

    // build linear system
    int nfixed = kp.size();

    Eigen::MatrixXd temp_mat = (Eigen::MatrixXd)laplacian_mat;
    Eigen::MatrixXd coeff_mat(laplacian_mat.rows()+nfixed, laplacian_mat.cols());
    Eigen::MatrixXd mat_b(laplacian_mat.rows()+nfixed, 3);
    for (int i = 0; i < temp_mat.rows(); i++) {
        coeff_mat.row(i) = temp_mat.row(i);
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j);
        }
    }

//    printf("Finish buding upper side of matrix!\n");

    for (int i = temp_mat.rows(); i<coeff_mat.rows(); i++) {
        int id = kp[i-temp_mat.rows()];
        coeff_mat.row(i).setZero();
        coeff_mat(i, id) = 1;
        mat_b(i, 0) = kvt[i-temp_mat.rows()].x();
        mat_b(i, 1) = kvt[i-temp_mat.rows()].y();
        mat_b(i, 2) = kvt[i-temp_mat.rows()].z();
    }

//    printf("Finish buding lower side of matrix!\n");

    Eigen::SparseMatrix<double> L_mat;
    L_mat.resize(coeff_mat.rows(), coeff_mat.cols());
    L_mat = coeff_mat.sparseView();

    Eigen::SparseMatrix<double> mat_a = L_mat.transpose() * L_mat;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(mat_a);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(solver.info() == Eigen::Success);
    if (solver.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

//    printf("Finish solver setup!\n");

//    start = clock();

    v = solver.solve(L_mat.transpose() * mat_b).eval();

//
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

//    printf("Fairing done!\n");

    if (partial_flag){
        int i = 0;
        for (auto it = vertex_candidates.begin(); it != vertex_candidates.end(); ++it) {
            raw_vertices[*it] = QVector3D(v(i, 0), v(i, 1), v(i, 2));
            ++i;
        }
    } else {
        for (int i=0; i<nvert; ++i) {
            raw_vertices[i] = QVector3D(v(i, 0), v(i, 1), v(i, 2));
        }
    }
//    printf("Copy done!\n");
    vertex_candidates.clear();

}


void LapDeform::localDeform(map<int, vec3>& shift_map, vector<int>& vertexIndices,
                     vector<vec3>& raw_vertices, vector<std::set<int>>& neighbors, vector<int>& fflag)
{
    // construct sub-mesh
    vector<vec3> sv, stv;
    vector<bool> bsfix, bvalid;
    vector<int> vid2sid, sid2vid, sf;

    bvalid.resize(raw_vertices.size(), false);
    vid2sid.resize(raw_vertices.size(), -1);

    int nf = vertexIndices.size()/3;
    for(int i=0; i<nf; i++)
    {
        int a, b, c;
        a = vertexIndices[i*3]-1;
        b = vertexIndices[i*3+1]-1;
        c = vertexIndices[i*3+2]-1;

        if(fflag[a]>=0&&fflag[b]>=0&&fflag[c]>=0)
        {
            sf.push_back(a); sf.push_back(b); sf.push_back(c);
            bvalid[a] = true; bvalid[b] = true; bvalid[c] = true;
        }
    }

    for(int i=0; i<raw_vertices.size(); i++)
    {
        if(bvalid[i])
        {
            sv.push_back(raw_vertices[i]);
            if(fflag[i]>0)
            {
                stv.push_back(raw_vertices[i]);
                bsfix.push_back(true);
            }
            else
            {
                stv.push_back(vec3(-1, -1, -1));
                bsfix.push_back(false);
            }
            vid2sid[i] = sv.size()-1;
            sid2vid.push_back(i);
        }
    }

    for(int i=0; i<sf.size(); i++)
    {
        sf[i] = vid2sid[sf[i]];
    }

    for(auto it = shift_map.begin(); it != shift_map.end(); ++it) {
        int id = vid2sid[it->first];
        if(id>=0)
        {
            stv[id] = it->second; bsfix[id] = true;
        }
    }

    cout<<sv.size()<<" "<<sf.size()<<endl;

    clock_t start, finish;
    double duration;
    start = clock();

    Eigen::MatrixXd v;
    Eigen::MatrixXd f;

    v.resize(sv.size(), 3);
    for (int i=0; i<sv.size(); ++i) {
        v(i, 0) = sv[i].x();
        v(i, 1) = sv[i].y();
        v(i, 2) = sv[i].z();
    }

    f.resize(sf.size()/3, 3);
    for (int i=0; i<sf.size()/3; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = sf[3*i+j];
        }
    }

    int nfixed = 0;
    vector<int> fixid;
    for(int i=0; i<bsfix.size(); i++)
    {
        if(bsfix[i])
        {
           fixid.push_back(i); nfixed++;
        }
    }

    cout<<"fixed "<<nfixed<<endl;

    Eigen::SparseMatrix<double> G, lapmat;
    igl::grad(v,f,G);
    Eigen::VectorXd dblA;
    igl::doublearea(v,f,dblA);
    const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
    lapmat = G.transpose() * T * G;

    Eigen::MatrixXd lap = lapmat*v;

    cout<<lapmat.rows()<<" "<<lapmat.cols()<<endl;

    Eigen::MatrixXd temp_mat = (Eigen::MatrixXd)lapmat;
    Eigen::MatrixXd coeff_mat(lapmat.rows()+nfixed, lapmat.cols());
    Eigen::MatrixXd mat_b(lapmat.rows()+nfixed, 3);
    for (int i = 0; i < temp_mat.rows(); i++) {
        coeff_mat.row(i) = temp_mat.row(i);
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j);
        }
    }

//    printf("Finish buding upper side of matrix!\n");

    float wt = 0.5;

    for (int i = temp_mat.rows(); i<coeff_mat.rows(); i++) {
        int id = fixid[i-temp_mat.rows()];
        coeff_mat.row(i).setZero();
        coeff_mat(i, id) = wt;
        mat_b(i, 0) = wt*stv[id].x();
        mat_b(i, 1) = wt*stv[id].y();
        mat_b(i, 2) = wt*stv[id].z();
    }

//    printf("Finish buding lower side of matrix!\n");

    Eigen::SparseMatrix<double> L_mat;
    L_mat.resize(coeff_mat.rows(), coeff_mat.cols());
    L_mat = coeff_mat.sparseView();

    Eigen::SparseMatrix<double> mat_a = L_mat.transpose() * L_mat;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(mat_a);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(solver.info() == Eigen::Success);
    if (solver.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

    v = solver.solve(L_mat.transpose() * mat_b).eval();

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i=0; i<sv.size(); i++) {
       raw_vertices[sid2vid[i]] = vec3(v(i, 0), v(i, 1), v(i, 2));
    }
    cout<<"han finished"<<endl;
}

void LapDeform::save_deform_info(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
                      std::vector<vec3>& handlePositions, vector<bool>& bfixed)
{
     ofstream handle(std::string(localfolder()+"\\models\\handleid.in").c_str());
     for(int i=0; i<handleIndices.size(); i++)
     {
         float wt = 1.0;
         if(bfixed[handleIndices[i]]) wt = 2.0;
         handle<<handleIndices[i]<<" "<<wt<<endl;
     }
     handle.close();

   /*  int nvert;
     int nface;
     Eigen::MatrixXd v, u, u_bc, v_bc;
     Eigen::MatrixXi f;

     nvert = verts.size();
     nface = vertIndices.size()/3;

     v.resize(nvert, 3);

     for (int i=0; i<nvert; ++i) {
         v(i, 0) = verts[i].x();
         v(i, 1) = verts[i].y();
         v(i, 2) = verts[i].z();
     }

     f.resize(nface, 3);

     for (int i=0; i<nface; i++) {
         for (int j=0; j<3; j++) {
             f(i, j) = vertIndices[3*i+j] - 1;
         }
     }

     Eigen::SparseMatrix<double> G, lapmat;
     igl::grad(v,f,G);
     Eigen::VectorXd dblA;
     igl::doublearea(v,f,dblA);
     const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
     lapmat = G.transpose() * T * G;

     ofstream lapm(std::string(localfolder()+"\\models\\lapDeform\\lapmat.tri").c_str());
     for (int k=0; k<lapmat.outerSize(); ++k)
       for (Eigen::SparseMatrix<double>::InnerIterator it(lapmat,k); it; ++it)
       {
         lapm<<it.row()<<" "<<it.col()<<" "<<it.value()<<endl;
       }
     lapm.close();*/
}

//void LapDeform::save_deform_info2(std::vector<vec3>& verts, std::vector<int>& vertIndices, std::vector<int>& handleIndices,
//                      std::vector<vec3>& handlePositions, vector<bool>& bfixed)
//{
//     ofstream handle(std::string(localfolder()+"\\models\\lapDeform\\handleid2.in").c_str());
//     for(int i=0; i<handleIndices.size(); i++)
//     {
//         float wt = 1.0;
//         if(bfixed[handleIndices[i]]) wt = 2.0;
//         handle<<handleIndices[i]<<" "<<wt<<endl;
//     }
//     handle.close();
//}

void LapDeform::handle_deform(vector<vec3>& verts, vector<int>& vertIndices, vector<int>& handleIndices,
                       vector<vec3>& handlePositions, vector<bool>& bfixed)
{
   /* ofstream handle(std::string(localfolder()+"\\models\\lapDeform\\handleid.in").c_str());
    for(int i=0; i<handleIndices.size(); i++)
    {
        float wt = 1.0;
        if(bfixed[handleIndices[i]]) wt = 2.0;
        handle<<handleIndices[i]<<" "<<wt<<endl;
    }
    handle.close();

    return;*/

    int nvert;
    int nface;
    Eigen::MatrixXd v, u, u_bc, v_bc;
    Eigen::MatrixXi f;

    nvert = verts.size();
    nface = vertIndices.size()/3;

    v.resize(nvert, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = verts[i].x();
        v(i, 1) = verts[i].y();
        v(i, 2) = verts[i].z();
    }

    f.resize(nface, 3);

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = vertIndices[3*i+j] - 1;
        }
    }

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::MatrixXd lap = lapmat*v;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    // build linear system
    int nfixed = handleIndices.size();

    cout<<"hanhanhan fixed"<<nfixed<<endl;

 //   Eigen::SparseMatrix<double,Eigen::RowMajor> coeff_mat(lapmat.rows()+nfixed, lapmat.cols());
  //  coeff_mat.reserve(Eigen::VectorXi::Constant(nvert+nfixed, 10));
    Eigen::MatrixXd mat_b(nvert+nfixed, 3);

 //   for(int k=0; k<lap_i.size(); k++)
 //   {
 //      double weight = 1.0;
 //      if(bmouth[lap_i[k]]) weight = 3.0;
 //      coeff_mat.insert(lap_i[k], lap_j[k]) = weight*lap_v[k];
 //   }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i = 0; i < lapmat.rows(); i++) {
        double weight = 1.0;
        if(blapfixed[i]) weight = 3.0;
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j)*weight;
        }
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

 //  double vweight = 1.5;
//    printf("Finish buding upper side of matrix!\n");

    for (int i = nvert; i<coffemat.rows(); i++) {
      //  int id = handleIndices[i-nvert];
     //   if(!bfixed[i]) vweight = 0.5f;
       // coeff_mat.row(i).setZero();
      //  coffemat.insert(i, id) = vweight;
        float vweight = handle_weight[i-nvert];
        mat_b(i, 0) = handlePositions[i-nvert].x()*vweight;
        mat_b(i, 1) = handlePositions[i-nvert].y()*vweight;
        mat_b(i, 2) = handlePositions[i-nvert].z()*vweight;
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(lapsolver.info() == Eigen::Success);
    if (lapsolver.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

    v = lapsolver.solve(Amat.transpose() * mat_b).eval();

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    printf("Fairing done!\n");

     for (int i=0; i<nvert; ++i) {
            verts[i] = vec3(v(i, 0), v(i, 1), v(i, 2));
        }
    cout<<v(0,0)<<" "<<v(100,1)<<" "<<v(200, 2)<<endl;
    printf("Copy done!\n");

}

void LapDeform::handle_deform_right(vector<vec3>& verts, vector<int>& vertIndices, vector<int>& handleIndices,
                       vector<vec3>& handlePositions, vector<bool>& bfixed)
{
   /* ofstream handle(std::string(localfolder()+"\\models\\lapDeform\\handleidright.in").c_str());
    for(int i=0; i<handleIndices.size(); i++)
    {
        float wt = 1.0;
        if(bfixed[handleIndices[i]]) wt = 2.0;
        handle<<handleIndices[i]<<" "<<wt<<endl;
    }
    handle.close();

    return;
*/
    int nvert;
    int nface;
    Eigen::MatrixXd v, u, u_bc, v_bc;
    Eigen::MatrixXi f;

    nvert = verts.size();
    nface = vertIndices.size()/3;

    v.resize(nvert, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = verts[i].x();
        v(i, 1) = verts[i].y();
        v(i, 2) = verts[i].z();
    }

    f.resize(nface, 3);

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = vertIndices[3*i+j] - 1;
        }
    }

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::MatrixXd lap = lapmatright*v;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    // build linear system
    int nfixed = handleIndices.size();

    cout<<"hanhanhan fixed"<<nfixed<<endl;

    Eigen::MatrixXd mat_b(nvert+nfixed, 3);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i = 0; i < lapmatright.rows(); i++) {
        double weight = 1.0;
        if(blapfixed[i]) weight = 3.0;
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j)*weight;
        }
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

 //  double vweight = 1.5;
//    printf("Finish buding upper side of matrix!\n");

    for (int i = nvert; i<coffematright.rows(); i++) {
        float vweight = handle_weight_right[i-nvert];
        mat_b(i, 0) = handlePositions[i-nvert].x()*vweight;
        mat_b(i, 1) = handlePositions[i-nvert].y()*vweight;
        mat_b(i, 2) = handlePositions[i-nvert].z()*vweight;
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(lapsolver_right.info() == Eigen::Success);
    if (lapsolver_right.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

    v = lapsolver_right.solve(AmatRight.transpose() * mat_b).eval();

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    printf("Fairing done!\n");

     for (int i=0; i<nvert; ++i) {
            verts[i] = vec3(v(i, 0), v(i, 1), v(i, 2));
        }
    cout<<v(0,0)<<" "<<v(100,1)<<" "<<v(200, 2)<<endl;
    printf("Copy done!\n");

}

void LapDeform::handle_deform_left(vector<vec3>& verts, vector<int>& vertIndices, vector<int>& handleIndices,
                       vector<vec3>& handlePositions, vector<bool>& bfixed)
{
  /*  ofstream handle(std::string(localfolder()+"\\models\\lapDeform\\handleidleft.in").c_str());
    for(int i=0; i<handleIndices.size(); i++)
    {
        float wt = 1.0;
        if(bfixed[handleIndices[i]]) wt = 2.0;
        handle<<handleIndices[i]<<" "<<wt<<endl;
    }
    handle.close();
    return;*/

    int nvert;
    int nface;
    Eigen::MatrixXd v, u, u_bc, v_bc;
    Eigen::MatrixXi f;

    nvert = verts.size();
    nface = vertIndices.size()/3;

    v.resize(nvert, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = verts[i].x();
        v(i, 1) = verts[i].y();
        v(i, 2) = verts[i].z();
    }

    f.resize(nface, 3);

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = vertIndices[3*i+j] - 1;
        }
    }

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::MatrixXd lap = lapmatleft*v;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    // build linear system
    int nfixed = handleIndices.size();

    cout<<"hanhanhan fixed"<<nfixed<<endl;


    Eigen::MatrixXd mat_b(nvert+nfixed, 3);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i = 0; i < lapmatleft.rows(); i++) {
        double weight = 1.0;
        if(blapfixed[i]) weight = 3.0;
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j)*weight;
        }
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

 //  double vweight = 1.5;
//    printf("Finish buding upper side of matrix!\n");

    for (int i = nvert; i<coffematleft.rows(); i++) {
        float vweight = handle_weight_left[i-nvert];
        mat_b(i, 0) = handlePositions[i-nvert].x()*vweight;
        mat_b(i, 1) = handlePositions[i-nvert].y()*vweight;
        mat_b(i, 2) = handlePositions[i-nvert].z()*vweight;
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(lapsolver_left.info() == Eigen::Success);
    if (lapsolver_left.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

    v = lapsolver_left.solve(AmatLeft.transpose() * mat_b).eval();

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    printf("Fairing done!\n");

     for (int i=0; i<nvert; ++i) {
            verts[i] = vec3(v(i, 0), v(i, 1), v(i, 2));
        }
    cout<<v(0,0)<<" "<<v(100,1)<<" "<<v(200, 2)<<endl;
    printf("Copy done!\n");

}


void LapDeform::handle_deform_all(vector<vec3>& verts, vector<int>& vertIndices, vector<int>& handleIndices,
                       vector<vec3>& handlePositions, vector<bool>& bfixed)
{
   /* ofstream handle(std::string(localfolder()+"\\models\\lapDeform\\handleid2.in").c_str());
    for(int i=0; i<handleIndices.size(); i++)
    {
        float wt = 1.0;
        if(bfixed[handleIndices[i]]) wt = 2.0;
        handle<<handleIndices[i]<<" "<<wt<<endl;
    }
    handle.close();*/


    int nvert;
    int nface;
    Eigen::MatrixXd v, u, u_bc, v_bc;
    Eigen::MatrixXi f;

    nvert = verts.size();
    nface = vertIndices.size()/3;

    v.resize(nvert, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = verts[i].x();
        v(i, 1) = verts[i].y();
        v(i, 2) = verts[i].z();
    }

    f.resize(nface, 3);

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = vertIndices[3*i+j] - 1;
        }
    }

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::MatrixXd lap = lapmat*v;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    // build linear system
    int nfixed = handleIndices.size();

    cout<<"hanhanhan fixed"<<nfixed<<endl;

 //   Eigen::SparseMatrix<double,Eigen::RowMajor> coeff_mat(lapmat.rows()+nfixed, lapmat.cols());
  //  coeff_mat.reserve(Eigen::VectorXi::Constant(nvert+nfixed, 10));
    Eigen::MatrixXd mat_b(nvert+nfixed, 3);

 //   for(int k=0; k<lap_i.size(); k++)
 //   {
 //      double weight = 1.0;
 //      if(bmouth[lap_i[k]]) weight = 3.0;
 //      coeff_mat.insert(lap_i[k], lap_j[k]) = weight*lap_v[k];
 //   }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i = 0; i < lapmat.rows(); i++) {
        double weight = 1.0;
        if(blapfixed[i]) weight = 3.0;
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j)*weight;
        }
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

 //  double vweight = 1.5;
//    printf("Finish buding upper side of matrix!\n");

    for (int i = nvert; i<coffemat.rows(); i++) {
      //  int id = handleIndices[i-nvert];
     //   if(!bfixed[i]) vweight = 0.5f;
       // coeff_mat.row(i).setZero();
      //  coffemat.insert(i, id) = vweight;
        float vweight = handle_weight[i-nvert];
        mat_b(i, 0) = handlePositions[i-nvert].x()*vweight;
        mat_b(i, 1) = handlePositions[i-nvert].y()*vweight;
        mat_b(i, 2) = handlePositions[i-nvert].z()*vweight;
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    assert(lapsolver.info() == Eigen::Success);
    if (lapsolver.info() != Eigen::Success) {
        printf("Solver failed!\n");
    }

    v = lapsolver.solve(Amat.transpose() * mat_b).eval();

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    printf("Fairing done!\n");

     for (int i=0; i<nvert; ++i) {
            verts[i] = vec3(v(i, 0), v(i, 1), v(i, 2));
        }
    cout<<v(0,0)<<" "<<v(100,1)<<" "<<v(200, 2)<<endl;
    printf("Copy done!\n");

}

void LapDeform::handle_deform_cuda(vector<vec3>& verts, vector<int>& vertIndices, vector<int>& handleIndices,
                       vector<vec3>& handlePositions)
{
    // --- Initialize cuSPARSE
    cusparseHandle_t handle;
    cusparseCreate(&handle);

    int nvert;
    int nface;
    Eigen::MatrixXd v, u, u_bc, v_bc;
    Eigen::MatrixXi f;

    nvert = verts.size();
    nface = vertIndices.size()/3;

    v.resize(nvert, 3);

    for (int i=0; i<nvert; ++i) {
        v(i, 0) = verts[i].x();
        v(i, 1) = verts[i].y();
        v(i, 2) = verts[i].z();
    }

    f.resize(nface, 3);

    for (int i=0; i<nface; i++) {
        for (int j=0; j<3; j++) {
            f(i, j) = vertIndices[3*i+j] - 1;
        }
    }

    clock_t start, finish;
    double duration;

    start = clock();

    Eigen::SparseMatrix<double> G, laplacian_mat;
    igl::grad(v,f,G);
    Eigen::VectorXd dblA;
    igl::doublearea(v,f,dblA);
    const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
    laplacian_mat = G.transpose() * T * G;

    Eigen::MatrixXd lap = laplacian_mat*v;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );


    // build linear system
    int nfixed = handleIndices.size();

    Eigen::SparseMatrix<double> L_mat;
    L_mat.resize(nvert+nfixed, nvert);

    laplacian_mat.makeCompressed();
    int* ptr = laplacian_mat.outerIndexPtr();
    int* col = laplacian_mat.innerIndexPtr();
    double* data = laplacian_mat.valuePtr();



    Eigen::MatrixXd mat_b(nvert+nfixed, 3);
    for (int i = 0; i < nvert; i++) {
        for (int j=0; j<3; j++) {
            mat_b(i, j) = lap(i, j);
        }
    }

    double weight = 1000;
//    printf("Finish buding upper side of matrix!\n");

    Eigen::SparseMatrix<double> mat_a = L_mat.transpose() * L_mat;



    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

/*
    Eigen::SparseMatrix<double,Eigen::RowMajor> L_mat(laplacian_mat);

    L_mat.makeCompressed();

    cusparseMatDescr_t descrA;	cusparseCreateMatDescr(&descrA);
    cusparseSetMatType		(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase	(descrA, CUSPARSE_INDEX_BASE_ZERO);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    int nfixed = handleIndices.size();

    int Nrows, nnz;
    Nrows = nvert; nnz = L_mat.nonZeros();

    std::cout<<nvert<<" "<<nnz<<std::endl;

    int* ptr = L_mat.outerIndexPtr();
    int* col = L_mat.innerIndexPtr();
    double* data = L_mat.valuePtr();

    std::cout<<nvert<<" "<<nnz<<std::endl;
    float weight = 1000;

    // --- Device side dense matrix
    double *d_A;			cudaMalloc(&d_A, (nnz+nfixed) * sizeof(*d_A));
    int *d_A_RowIndices;	cudaMalloc(&d_A_RowIndices, (Nrows +nfixed + 1) * sizeof(*d_A_RowIndices));
    int *d_A_ColIndices;    cudaMalloc(&d_A_ColIndices, (nnz+nfixed) * sizeof(*d_A_ColIndices));

    std::cout<<nfixed<<" "<<nnz<<std::endl;

    for(int i=0; i<nnz; i++)
    {
        d_A[i] = data[i]; d_A_ColIndices[i] = col[i];
        std::cout<<d_A[i]<<" i "<<d_A_ColIndices[i]<<std::endl;
    }
    for(int i=0; i<Nrows; i++)
    {
        d_A_RowIndices[i] = ptr[i];
        std::cout<<" i "<<d_A_RowIndices[i]<<std::endl;
    }

    for(int i=0; i<nfixed; i++)
    {
        d_A[nnz+i] = weight;
        d_A_ColIndices[nnz+i] = handleIndices[i];
        d_A_RowIndices[Nrows+i] = nnz + i;
    }

    d_A_RowIndices[Nrows+nfixed] = nnz + nfixed;

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    // --- Host side dense matrix
    double *h_A = (double *)malloc((nnz+nfixed) * sizeof(*h_A));
    int *h_A_RowIndices = (int *)malloc( (Nrows +nfixed + 1) * sizeof(*h_A_RowIndices));
    int *h_A_ColIndices = (int *)malloc((nnz+nfixed) * sizeof(*h_A_ColIndices));



    cudaMemcpy(h_A, d_A, (nnz+nfixed)*sizeof(*h_A), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_A_RowIndices, d_A_RowIndices, (Nrows +nfixed + 1) * sizeof(*h_A_RowIndices), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_A_ColIndices, d_A_ColIndices, (nnz+nfixed) * sizeof(*h_A_ColIndices), cudaMemcpyDeviceToHost);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );

    for (int i = 0; i < nnz; ++i) printf("A[%i] = %.0f ", i, h_A[i]); printf("\n");

    for (int i = 0; i < (Nrows + 1); ++i) printf("h_A_RowIndices[%i] = %i \n", i, h_A_RowIndices[i]); printf("\n");

    for (int i = 0; i < nnz; ++i) printf("h_A_ColIndices[%i] = %i \n", i, h_A_ColIndices[i]);


    // --- CUDA solver initialization
    cusolverSpHandle_t solver_handle;
    cusolverSpCreate(&solver_handle);

    for(int k=0; k<1; k++)
    {
        // --- Allocating and defining dense host and device data vectors
        double *h_y	= (double *)malloc((Nrows+nfixed) * sizeof(double));
        for(int i=0; i<Nrows; i++)
        {
            h_y[i] = lap(i, k);
        }
        for(int i=0; i<nfixed; i++)
        {
            h_y[i+Nrows] = weight*handlePositions[i][k];
        }

        double *d_y;		cudaMalloc(&d_y, (Nrows+nfixed) * sizeof(double));
        cudaMemcpy(d_y, h_y, (Nrows+nfixed) * sizeof(double), cudaMemcpyHostToDevice);

        // --- Allocating the host and device side result vector
        double *h_x	= (double *)malloc(Nrows * sizeof(double));
        double *d_x;		cudaMalloc(&d_x, Nrows * sizeof(double));


        int rankA;
        int *p = (int *)malloc(Nrows * sizeof(int));
        double min_norm;

        cusolverSpDcsrlsqvqrHost(solver_handle, Nrows+nfixed, Nrows, nnz, descrA, h_A, h_A_RowIndices, h_A_ColIndices, h_y, 0.000001, &rankA, h_x, p, &min_norm);

        printf("Showing the results...\n");
        for (int i = 0; i < Nrows; i++) printf("%f\n", h_x[i]);
    }*/

}


void LapDeform::lapDeform_cuda(const std::map<int, QVector3D> & shift_map,
                     const std::vector<int> & vertexIndices,
                     std::vector<QVector3D> & raw_vertices,
                     const std::vector<std::set<int> > & neighbors,
                     bool debug_flag,
                     bool partial_flag,
                     bool shift_flag){



}
