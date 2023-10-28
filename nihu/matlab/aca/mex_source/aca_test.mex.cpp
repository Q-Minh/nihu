#include <complex>
#include <Eigen/Core>
#include <mex.h>
#include "nihu/util/mex_matrix.hpp"
#include "nihu/aca/aca.hpp"

class HelmholtzMatrix
{
public:
	typedef std::complex<double> Scalar;

	HelmholtzMatrix(
		Eigen::Matrix<double, Eigen::Dynamic, 2> const &xs,
		Eigen::Matrix<double, Eigen::Dynamic, 2> const &xr,
		std::complex<double> const &k) : m_xs(xs), m_xr(xr), m_k(k)
	{
	}

	Scalar operator()(int i, int j) const
	{
		double r = (m_xr.row(i) - m_xs.row(j)).norm();
		return exp(std::complex<double>(0.,-1.)*m_k*r)/r;
	}

private:
	Eigen::Matrix<double, Eigen::Dynamic, 2> m_xs;
	Eigen::Matrix<double, Eigen::Dynamic, 2> m_xr;
	std::complex<double> m_k;
};


// [resp, outranks] = aca_test(k, xs, xr, CS, CR, B_far, eps, R, exc)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	std::complex<double> k(*mxGetPr(prhs[0]), *mxGetPi(prhs[0]));
	mex::real_matrix<double> xs(prhs[1]), xr(prhs[2]);
	mex::real_matrix<int> CS(prhs[3]), CR(prhs[4]), B_far(prhs[5]);
	double eps = mxGetScalar(prhs[6]);
	int R = *((int *)mxGetData(prhs[7]));
	mex::complex_matrix<double> exc(prhs[8]);

	mex::complex_matrix<double> res(xr.rows(), 1, plhs[0]);
	for (int i = 0; i < res.rows(); ++i) res(i,0) = 0.0;
	mex::real_matrix<int> outRanks(B_far.rows(), 1, plhs[1]);

	HelmholtzMatrix HM(xs, xr, k);

	ACA aca;
	aca.multiply(HM, CR, CS, B_far, exc, res, eps, R, outRanks);
}

