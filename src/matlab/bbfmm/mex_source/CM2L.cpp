#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "mex.h"

typedef double Real;
typedef int Index;
typedef Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> > MapMatrix;
typedef Eigen::Map<Eigen::Matrix<Index, Eigen::Dynamic, 1> > MapVector;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
// Local = M2L(source, receiver, receiveridx, U, V, ranks, multi)
{
    mxArray const *source = prhs[0];
    mxArray const *receiver = prhs[1];
    mxArray const *receiveridx = prhs[2];
    mxArray const *M2L_U = prhs[3];
    mxArray const *M2L_V = prhs[4];
    mxArray const *M2L_ranks = prhs[5];
    mxArray const *multi = prhs[6];
    
    int N = mxGetM(M2L_U);
    int maxR = mxGetN(M2L_U)/343;
    int nClusters = mxGetN(multi);
    int K = mxGetM(source);
    
    plhs[0] = mxCreateNumericMatrix(N, nClusters, mxDOUBLE_CLASS, mxREAL);
    mxArray *local = plhs[0];
    
    MapVector s((int *)mxGetData(source), K);
    MapVector r((int *)mxGetData(receiver), K);
    MapVector ridx((int *)mxGetData(receiveridx), K);
    
    MapMatrix m(mxGetPr(multi), N, nClusters);
    MapMatrix l(mxGetPr(local), N, nClusters);
    std::vector<MapMatrix> U, V;
    U.reserve(343);
    V.reserve(343);
    for (int i = 0; i < 343; ++i)
    {
        U.push_back(MapMatrix(mxGetPr(M2L_U)+(N*maxR*i), N, maxR));
        V.push_back(MapMatrix(mxGetPr(M2L_V)+(N*maxR*i), N, maxR));
    }
    MapVector ranks((int *)mxGetData(M2L_ranks), 343);
    
    // the actual operation
    l.setZero();
    for (int i = 0; i < s.rows(); ++i)
    {
        int rx = ridx(i);
        int rnk = ranks(rx);
        l.col(r(i)) += U[rx].leftCols(rnk) *
                ( V[rx].leftCols(rnk).adjoint() * m.col(s(i)) );
    }
}

