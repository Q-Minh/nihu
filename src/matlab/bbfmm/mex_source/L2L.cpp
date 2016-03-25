#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "mex.h"

typedef double Real;
typedef int Index;
typedef Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
typedef Eigen::Map<Eigen::Matrix<Index, Eigen::Dynamic, 1> > MapVector;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
// chLocal = L2L(local, source, receiver, receiveridx, L2L, size)
{
    mxArray const *local = prhs[0];
    mxArray const *source= prhs[1];
    mxArray const *receiver = prhs[2];
    mxArray const *receiveridx = prhs[3];
    mxArray const *L2L = prhs[4];
    mxArray const *size = prhs[5];
    
    int N = mxGetM(local);
    int nClusters = mxGetN(local);
    int nChildren = *(int *)mxGetData(size);
    int K = mxGetM(source);
    
    plhs[0] = mxCreateNumericMatrix(N, nChildren, mxDOUBLE_CLASS, mxREAL);
    mxArray *childLocal = plhs[0];
    
    MapVector s((int *)mxGetData(source), K);
    MapVector r((int *)mxGetData(receiver), K);
    MapVector ridx((int *)mxGetData(receiveridx), K);
    
    MapMatrix l(mxGetPr(local), N, nClusters);
    MapMatrix cl(mxGetPr(childLocal), N, nChildren);

    std::vector<MapMatrix> l2l;
    l2l.reserve(8);
    for (int i = 0; i < 8; ++i)
        l2l.push_back(MapMatrix(mxGetPr(L2L)+(N*N*i), N, N));
    
    // the actual operation
    cl.setZero();
    for (int i = 0; i < K; ++i)
        cl.col(r(i)) += l2l[ridx(i)] * l.col(s(i));
}
