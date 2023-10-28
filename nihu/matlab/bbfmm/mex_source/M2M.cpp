#include <Eigen/Core>
#include <iostream>
#include <vector>
#include "mex.h"

typedef double Real;
typedef int Index;
typedef Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
typedef Eigen::Map<Eigen::Matrix<Index, Eigen::Dynamic, 1> > MapVector;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
// Local = M2L(source, receiver, receiveridx, M2L, multi)
{
    mxArray const *multi = prhs[0];
    mxArray const *father= prhs[1];
    mxArray const *fatheridx = prhs[2];
    mxArray const *M2M = prhs[3];
    mxArray const *size = prhs[4];
    
    int N = mxGetM(multi);
    int nClusters = mxGetN(multi);
    int nParent = *(int *)mxGetData(size);
    
    plhs[0] = mxCreateNumericMatrix(N, nParent, mxDOUBLE_CLASS, mxREAL);
    mxArray *parentMulti = plhs[0];
    
    MapVector f((int *)mxGetData(father), nClusters);
    MapVector fidx((int *)mxGetData(fatheridx), nClusters);
    
    MapMatrix m(mxGetPr(multi), N, nClusters);
    MapMatrix pm(mxGetPr(parentMulti), N, nParent);

    std::vector<MapMatrix> m2m;
    m2m.reserve(8);
    for (int i = 0; i < 8; ++i)
        m2m.push_back(MapMatrix(mxGetPr(M2M)+(N*N*i), N, N));
    
    // the actual operation
    pm.setZero();
    for (int i = 0; i < nClusters; ++i)
        pm.col(f(i)) += m2m[fidx(i)] * m.col(i);
}
