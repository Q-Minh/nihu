#include <Eigen/Dense>
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
	mxArray const *source = prhs[0];
	mxArray const *receiver = prhs[1];
	mxArray const *receiveridx = prhs[2];
	mxArray const *M2L = prhs[3];
	mxArray const *multi = prhs[4];
	
	int N = mxGetM(multi);
	int nClusters = mxGetN(multi);
	int K = mxGetM(source);
	
	plhs[0] = mxCreateNumericMatrix(N, nClusters, mxDOUBLE_CLASS, mxREAL);
	mxArray *local = plhs[0];
	
	MapVector s((int *)mxGetData(source), K);
	MapVector r((int *)mxGetData(receiver), K);
	MapVector ridx((int *)mxGetData(receiveridx), K);
	
	MapMatrix m(mxGetPr(multi), N, nClusters);
	MapMatrix l(mxGetPr(local), N, nClusters);
	std::vector<MapMatrix> m2l;
    m2l.reserve(343);
	for (int i = 0; i < 343; ++i)
		m2l.push_back(MapMatrix(mxGetPr(M2L)+(N*N*i), N, N));
		
	// the actual operation
	l.setZero();
	for (int i = 0; i < s.rows(); ++i)
		l.col(r(i)) += m2l[ridx(i)] * m.col(s(i));
}

