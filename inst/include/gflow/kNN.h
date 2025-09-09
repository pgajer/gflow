#ifndef MSC2_KNN_H_
#define MSC2_KNN_H_

#include <ANN/ANN.h>

struct knn_result_t {
    int* indices;
    double* distances;
};

SEXP kNN_search(ANNpointArray data_pts, int nrX, int ncX, int k);

#endif // MSC2_KNN_H_
