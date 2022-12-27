//
// Created by 杜扬 on 2022/4/24.
//

#ifndef CSPIKE_SPIKE_SKETCH_H
#define CSPIKE_SPIKE_SKETCH_H

#include "utils/MurmurHash3.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

class spike_sketch {
public:
    spike_sketch(int n, int p, int ncodes, uint32_t seed);
    int n;//cell数量
    int p;//更新概率
    int sign;
    double q;//
    int ncodes;
    uint32_t seed;
    uint32_t seed2;
    int* S;
    int* E;
//    bool* T;
    int vlimit;
    int slimit;
    void updateT();
    void update(int key);
    double query(double alpha0,double alpha1,double beta0,double beta1,double coe);
    bool valid();
    bool tension(int j);
    void adjust(int j, int v);
    void extended_inc(int j, int v);
    int rho(uint64_t x);

    double getF_R_j(int j,double alpha0,double alpha1,double beta0,double beta1);
    double getZ_k_0(double k,double alpha0);
    double getZ_k_1(double k,double alpha1);
    double getZ_k(double k);
    double getZ_k_minus(double k,double alpha0,double alpha1);

    int getNumOfStages();
    int getStageWidth();
};


#endif //CSPIKE_SPIKE_SKETCH_H
