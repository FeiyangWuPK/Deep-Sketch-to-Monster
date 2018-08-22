#ifndef CAFFE_PARSER_H
#define CAFFE_PARSER_H

#include "caffeHeader.h"
#include "caffe/caffe.hpp"
#include "caffe/util/io.hpp"
#include "caffe/blob.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/data_transformer.hpp"
#include <string>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include "folder.h"

#define TENSOR_SIZE 34530
#define U_SIZE 50
#define V_SIZE 16

using namespace caffe;

Net<double> gesture_net(localfolder()+"\\models\\caffe\\gesture.prototxt", TEST);
Net<double> net_u(localfolder()+"\\models\\caffe\\regression_u50.prototxt", TEST);
Net<double> net_v(localfolder()+"\\models\\caffe\\regression_v16.prototxt", TEST);

// ################ Gesture Net ################
int gesture_run(const std::string & imagefile){
    //get datum
    Datum datum;
    if (!ReadImageToDatum(imagefile, 1, 227, 227, false, &datum)) { //is_color = false
        LOG(ERROR) << "Error during file reading";
        return gesture_run(imagefile);
    }

    //set image transformer
    TransformationParameter transformationParameter;
    transformationParameter.set_mirror(false);
    transformationParameter.set_crop_size(227);
    transformationParameter.add_mean_value(200);
    DataTransformer<double> dataTransformer(transformationParameter, TEST);

    //get the blob
    Blob<double>* blob = new Blob<double>(1, datum.channels(), datum.height(), datum.width());

    //transform image
    dataTransformer.Transform(datum, blob);

    //fill the vector
    vector<Blob<double>*> bottom;
    bottom.push_back(blob);
    double loss_type = 0.0;


    //forward the net
    const vector<Blob<double>*>& netbottom = gesture_net.Forward(bottom, &loss_type);

    float max_result = -1;
    int max_pos = -1;

    for (int i = 0; i < 14; ++i){
        float current_max = netbottom[0]->cpu_data()[i];
        if (current_max > max_result){
            max_result = current_max;
            max_pos = i;
        }
    }

    if (max_result < 0.2) return 15;
    std::cout << max_pos << std::endl;
    return max_pos + 1;
}


// ################ UV Net ################
Blob<double>* transblob;
Blob<double> middleblob, resultblob;


void caffe_setup(){
    //Caffe::set_mode(Caffe::CPU);
    Caffe::set_mode(Caffe::GPU);
    Caffe::SetDevice(0);

    // load caffemodel
    net_u.CopyTrainedLayersFrom(localfolder()+"\\models\\caffe\\regression_u50.caffemodel");
    net_v.CopyTrainedLayersFrom(localfolder()+"\\models\\caffe\\regression_v16.caffemodel");

    //load transform matrix into Blob<double>* transblob
    double* C = (double *)malloc(TENSOR_SIZE*U_SIZE*V_SIZE*sizeof(double));
    FILE * pFile = fopen (std::string(localfolder()+"\\models\\caffe\\model_tensor.bin").c_str(), "rb");
    fread (C , sizeof(double), TENSOR_SIZE*U_SIZE*V_SIZE, pFile);
    transblob = new Blob<double>(1, TENSOR_SIZE*U_SIZE*V_SIZE, 1, 1);

    BlobProto transblob_proto;
    transblob_proto.set_num(1);
    transblob_proto.set_channels(TENSOR_SIZE*U_SIZE*V_SIZE);
    transblob_proto.set_height(1);
    transblob_proto.set_width(1);
    for (int i = 0; i < TENSOR_SIZE*U_SIZE*V_SIZE; ++i) {
        transblob_proto.add_data(0.);
        transblob_proto.set_data(i, transblob_proto.data(i) + C[i]);
    }
    free(C);
    transblob->FromProto(transblob_proto);
    std::vector<int> transhape = {1, TENSOR_SIZE*U_SIZE*V_SIZE};
    transblob->Reshape(transhape);

    //set buffer blob
    std::vector<int> middleblobshape = {1, TENSOR_SIZE*U_SIZE};
    std::vector<int> resultshape = {1, TENSOR_SIZE};
    middleblob.Reshape(middleblobshape);
    resultblob.Reshape(resultshape);

    gesture_net.CopyTrainedLayersFrom(localfolder()+"\\models\\caffe\\gesture.caffemodel");
}


const double* caffe_run(const std::string & imagefile){

    clock_t begin_time = clock();

    //get datum
    Datum datum;
    if (!ReadImageToDatum(imagefile, 1, 256, 256, false, &datum)) { //is_color = false
        LOG(ERROR) << "Error during file reading";
    }

    printf("    Image Reading Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


    //set image transformer
    TransformationParameter transformationParameter;
    transformationParameter.set_mirror(false);
    transformationParameter.set_crop_size(227);
    transformationParameter.add_mean_value(200);
    DataTransformer<double> dataTransformer(transformationParameter, TEST);

    //get the blob
    Blob<double>* blob = new Blob<double>(1, datum.channels(), 227, 227);

    //transform image
    dataTransformer.Transform(datum, blob);

    //fill the vector
    vector<Blob<double>*> bottom;
    bottom.push_back(blob);
    double loss_type = 0.0;

    printf("    Image Layer Setup Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


    //forward the net
    const vector<Blob<double>*>& netbottom_u = net_u.Forward(bottom, &loss_type);
    const vector<Blob<double>*>& netbottom_v = net_v.Forward(bottom, &loss_type);

    printf("    Net Forward Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


//    for (int t = 0; t < U_SIZE; ++t){
//        printf("u_%d: %f\n", t, netbottom_u[0]->cpu_data()[t]);
//    }
//    for (int t = 0; t < V_SIZE; ++t){
//        printf("v_%d: %f\n", t, netbottom_v[0]->cpu_data()[t]);
//    }

    // start calculating
    caffe_gpu_gemv(
        CblasNoTrans,
        TENSOR_SIZE*U_SIZE, V_SIZE,
        1.0,
        transblob->gpu_data(),
        netbottom_v[0]->gpu_data(),
        0.0, middleblob.mutable_gpu_data());

    caffe_gpu_gemv(
        CblasNoTrans,
        TENSOR_SIZE, U_SIZE,
        1.0,
        middleblob.gpu_data(),
        netbottom_u[0]->gpu_data(),
        0.0, resultblob.mutable_gpu_data());

    printf("    New Vectors Calculating Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();

    return resultblob.cpu_data();
}

const double* caffe_view2_run(const std::string & imagefile1, const std::string & imagefile2){

    clock_t begin_time = clock();

    //get datum
    Datum datum1;
    if (!ReadImageToDatum(imagefile1, 1, 256, 256, false, &datum1)) { //is_color = false
        LOG(ERROR) << "Error during file reading";
    }

    Datum datum2;
    if (!ReadImageToDatum(imagefile2, 1, 256, 256, false, &datum2)) { //is_color = false
        LOG(ERROR) << "Error during file reading";
    }

    printf("    Image Reading Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


    //set image transformer
    TransformationParameter transformationParameter;
    transformationParameter.set_mirror(false);
    transformationParameter.set_crop_size(227);
    transformationParameter.add_mean_value(200);
    DataTransformer<double> dataTransformer(transformationParameter, TEST);

    printf("  blob1 channels %d\n", datum1.channels());
    printf("  blob1 channels %d\n", datum2.channels());

    //get the blob
    Blob<double>* blob1 = new Blob<double>(1, datum1.channels(), 227, 227);
    Blob<double>* blob2 = new Blob<double>(1, datum2.channels(), 227, 227);

    //transform image
    dataTransformer.Transform(datum1, blob1);
    dataTransformer.Transform(datum2, blob2);

    printf(" data transform %d\n", datum1.height());
    printf(" data transform %d\n", datum2.width());

    //fill the vector
    vector<Blob<double>*> bottom1;
    bottom1.push_back(blob1);
    bottom1.push_back(blob2);

    vector<Blob<double>*> bottom2;
    bottom2.push_back(blob1);
    double loss_type = 0.0;

    printf("    Image Layer Setup Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


    //forward the net
    const vector<Blob<double>*>& netbottom_u = net_u.Forward(bottom1, &loss_type);
    const vector<Blob<double>*>& netbottom_v = net_v.Forward(bottom2, &loss_type);

    printf("    Net Forward Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();


   /* for (int t = 0; t < U_SIZE; ++t){
        printf("u_%d: %f\n", t, netbottom_u[0]->cpu_data()[t]);
    }
    for (int t = 0; t < V_SIZE; ++t){
        printf("v_%d: %f\n", t, netbottom_v[0]->cpu_data()[t]);
    }
*/
    // start calculating
    caffe_gpu_gemv(
        CblasNoTrans,
        TENSOR_SIZE*U_SIZE, V_SIZE,
        1.0,
        transblob->gpu_data(),
        netbottom_v[0]->gpu_data(),
        0.0, middleblob.mutable_gpu_data());

    caffe_gpu_gemv(
        CblasNoTrans,
        TENSOR_SIZE, U_SIZE,
        1.0,
        middleblob.gpu_data(),
        netbottom_u[0]->gpu_data(),
        0.0, resultblob.mutable_gpu_data());

    printf("    New Vectors Calculating Time: %f\n", float( clock () - begin_time ) /  CLOCKS_PER_SEC); begin_time = clock ();

    return resultblob.cpu_data();
}

#endif
