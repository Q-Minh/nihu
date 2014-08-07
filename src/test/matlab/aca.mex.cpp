#include "util/mex_matrix.hpp"

#include "aca/aca.hpp"


template <class LocMatrix>
class Laplace
{
public:
	typedef double Scalar;

	Laplace(LocMatrix const &xr, LocMatrix const &xs) : xr(xr), xs(xs)
	{
	}

	double operator()(int i, int j) const
	{
		return 1./(xr.row(i) - xs.row(j)).norm();
	}

private:
	LocMatrix const &xr;
	LocMatrix const &xs;
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	mex::real_matrix<double> xs(prhs[0]), xr(prhs[1]);
	mex::real_matrix<int> source_clusters(prhs[2]), receiver_clusters(prhs[3]), blocks(prhs[4]);
	mex::real_matrix<double> input(prhs[5]);

	mex::real_matrix<double> output(input.rows(), 1, plhs[0]);
	output.setZero();

	Laplace<mex::real_matrix<double> > M(xr, xs);

	ACA aca;
	aca.multiply(M, receiver_clusters, source_clusters, blocks, input, output, 1e-3, 50);
}

