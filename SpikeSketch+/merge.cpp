#include "spike_sketch_extend.h"
using namespace std;
void max_calcu(int n,int &max,int &sub_max){
    if(n>max){
        sub_max=max;
        max=n;
    }
    else if(n>sub_max){
        sub_max=n;
    }
}
void merge(vector<vector<spike_sketch>> spikeSketchArray,vector<spike_sketch>& merged_ss){
    int numOfBuckets = spikeSketchArray[0].size();
    int numOfCells=spikeSketchArray[0][0].n;
    int numOfSS=spikeSketchArray.size();
//    vector<int[2]> A(numOfCells);
    for(int i=0;i<numOfBuckets;i++){
        for(int j=0;j<numOfCells;j++){
            int max=0;
            int sub_max=0;
            for(int k=0;k<numOfSS;k++){
                int S=spikeSketchArray[k][i].S[j];
                int E=spikeSketchArray[k][i].E[j];
                if(spikeSketchArray[k][i].tension(j) && E==1){
                    max_calcu(S,max,sub_max);
                    max_calcu(S-1,max,sub_max);
                }
                else if(!spikeSketchArray[k][i].tension(j) && E==0){
                    continue;
                }
                else{
                    max_calcu(S,max,sub_max);
                }
            }
            merged_ss[i].adjust(j,max);
            merged_ss[i].adjust(j,sub_max);
        }
    }
}
