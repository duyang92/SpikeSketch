#ifndef _SPIKEUTIL_H_
#define _SPIKEUTIL_H_

#include <iostream>
#include "MurmurHash3.h"
#include <list>
#include <numeric>

using namespace std;

template <typename T>
double stats_mean(T data[], int n)
{
    if (n <= 0)
    {
        std::cerr << "wrong parameter!" << std::endl;
        return -1;
    }

    double mean = 0.0;
    for (int i = 0; i < n; i++)
        mean += data[i];
    mean /= n;

    return mean;
}

//template <typename T>
//double stats_variance(T data[], int n)
//{
//    if (n <= 0)
//    {
//        std::cerr << "wrong parameter!" << std::endl;
//        return -1;
//    }
//
//    double mean = 0.0;
//    for (int i = 0; i < n; i++)
//        mean += data[i];
//    mean /= n;
//
//    double variance = 0.0;
//    for (int i = 0; i < n; i++)
//        variance += pow(data[i] - mean, 2);
//    variance /= n;
//
//    return variance;
//}

template <typename T>
double stats_variance(T data[],double mean, int n)
{
    if (n <= 0)
    {
        std::cerr << "wrong parameter!" << std::endl;
        return -1;
    }

    double variance = 0.0;
    for (int i = 0; i < n; i++)
        variance += pow(data[i] - mean, 2);
    variance /= n;

    return variance;
}
#endif