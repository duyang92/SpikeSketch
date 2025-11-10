#include <iostream>
#include "utils/MurmurHash3.h"
#include "spike_sketch_extend.h"

#include <fstream>
#include <ctime>
#include <random>
#include <iostream>
#include <sstream>

using namespace std;

static mt19937 rng(time(NULL));

double mutiBktQuery(vector<spike_sketch> &spikeSketchArray, double alpha0, double alpha1, double beta0, double beta1,
                    double coe) {
    double threshold = 0.04;

    int numOfBuckets = spikeSketchArray.size();

    int numOfStages = 0;
    int stageWidth = 0;
    if (numOfBuckets > 0) {
        numOfStages = spikeSketchArray[0].getNumOfStages();
        stageWidth = spikeSketchArray[0].getStageWidth();
    }

    vector<double> emptyRatios;
    vector<int> numOfRegsInStages;
    int fzSum = 0;

    for (int stageIdx = 0; stageIdx < numOfStages; stageIdx++) {
        double totalRegs = 0;
        double unknownRegs[32] = {0};
        double rankRegs[32] = {0};
        double rankNonEmptyRegs[32] = {0};

        for (int bucketIdx = 0; bucketIdx < numOfBuckets; bucketIdx++) {
            spike_sketch ss = spikeSketchArray[bucketIdx];

            for (int cellIdx = stageIdx * stageWidth;
                 cellIdx < (stageIdx + 1) * stageWidth && cellIdx < ss.n - 1; cellIdx++) {
                totalRegs++;
                if (!ss.tension(cellIdx) || ss.E[cellIdx] == 1) {
                    rankRegs[ss.S[cellIdx]] += 1;
                    if (ss.S[cellIdx] != 0) {
                        rankNonEmptyRegs[ss.S[cellIdx]] += 1;
                    }
                } else {
                    if (ss.S[cellIdx] != 1) {
                        unknownRegs[ss.S[cellIdx]] += 1;
                    } else {
                        rankRegs[ss.S[cellIdx]] += 1;
                    }
                }
            }
        }

        double smallerRegs = rankRegs[0] + rankRegs[1];
        double smallerNonEmpty = rankNonEmptyRegs[0] + rankNonEmptyRegs[1];
        double smallerNonEmptyRatio;
        for (int r = 2; r < 32; r++) {
            smallerNonEmptyRatio = smallerNonEmpty / smallerRegs;
            smallerRegs += unknownRegs[r];
            smallerRegs += rankRegs[r];
            smallerNonEmpty += smallerNonEmptyRatio * unknownRegs[r];
            smallerNonEmpty += rankNonEmptyRegs[r];
        }

        double emptyRatio = 1 - (double) smallerNonEmpty / smallerRegs;

        if (emptyRatio >= threshold) {
            fzSum++;
            emptyRatios.push_back(emptyRatio);
        }
    }


    if (fzSum > 2) {//Estimate for small flow
        double sum = 0;
        for (int stageIdx = 0; stageIdx < numOfStages; stageIdx++) {
            if (emptyRatios[stageIdx] >= threshold) {
                double o_z = double(numOfStages - 1 - stageIdx) / numOfStages;
                sum += -1.0 * pow(4, o_z) * numOfBuckets * 5 * log(emptyRatios[stageIdx]);
            }
        }
        double result = ((double) numOfStages / fzSum) * sum;
        return result;
    } else {
        double result = 0;
        for (int bucketIdx = 0; bucketIdx < numOfBuckets; bucketIdx++) {
            spike_sketch ss = spikeSketchArray[bucketIdx];
            result += ss.query(alpha0, alpha1, beta0, beta1, coe);
        }
        return result;
    }
}

double mutiBktQuery1(vector<spike_sketch> &spikeSketchArray, double alpha0, double alpha1, double beta0, double beta1,
                    double coe) {

    int numOfBuckets = spikeSketchArray.size();

    double result = 0;
    for (int bucketIdx = 0; bucketIdx < numOfBuckets; bucketIdx++) {
        result += spikeSketchArray[bucketIdx].query(alpha0, alpha1, beta0, beta1, coe);
    }
    return result;

}


int main() {
    srand(time(NULL));
    
    int nf = 50;//Actual Cardinality

    uint64_t key;
    int n = 20;//Number of cells in a single bucket
    int ncode = 4;//Number of bits in a single cell
    int p = 12;

    double alpha0 = 0.1;
    double alpha1 = 0.88;
    double beta0 = 1.12;
    double beta1 = 1.46;
    double myCoe = 0.573;//Coefficient of correction

    int memory = 8*64;
    int numOfBuckets = memory>>6;
    
    time_t seed = rand() + 1316561;

    vector<spike_sketch> spikeSketchArray;
    for (int bktIdx = 0; bktIdx < numOfBuckets; bktIdx++) {
        spikeSketchArray.push_back(spike_sketch(n, p, ncode, seed));
    }
    for (key = 1; key <= nf; key++) {
        uint32_t tempInt32 = 0;
        MurmurHash3_x86_32(&key, 4, seed + 231321, &tempInt32);
        spikeSketchArray[tempInt32 % numOfBuckets].update(key);

    double result = 0.0;
    for (int bucketIdx = 0; bucketIdx < numOfBuckets; bucketIdx++) {
        if (!spikeSketchArray[bucketIdx].getCompressed()) {
            result += spikeSketchArray[bucketIdx].query(alpha0, alpha1, beta0, beta1, myCoe);
            spikeSketchArray[bucketIdx] = spike_sketch(n, p, ncode, seed);
            spikeSketchArray[bucketIdx].setCompressed();
        }
    }
    result += mutiBktQuery(spikeSketchArray, alpha0, alpha1, beta0, beta1, myCoe);
            
    char traceFilePath[100];
    ofstream logFile;
    sprintf(traceFilePath, "result.txt");
    logFile.open(traceFilePath, ios::out);
    logFile<<"Input:\n\tActuall Cardinality:"<<nf<<"\n\tNumber of Buckets: "<<numOfBuckets<<endl;
    logFile << "######################################" << endl;
    logFile<<"Parameters of SpikeSketch:\n\t";
    logFile << "ncodes: " << ncode << ", ";
    logFile << "p: " << p << ", ";
    logFile << "ncell: " << n << endl;
    logFile << "\talpha0: " << alpha0 << ", ";
    logFile << "alpha1: " << alpha1 << endl;
    logFile << "\tbeta0: " << beta0 << ", ";
    logFile << "beta1: " << beta1 << endl;
    logFile << "\tcoe: " << myCoe << endl;
    logFile << "######################################" << endl;
    logFile << "Output:\n\t" << "Estimate Cardinality: " << result << "\n\tMemory Usage: " << memory;
    logFile.close();
    return 0;

}

