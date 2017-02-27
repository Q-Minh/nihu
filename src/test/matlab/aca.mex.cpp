#include "util/mex_matrix.hpp"

#include "aca/aca.hpp"


/** \brief function computing the Newton kernel matrix
 * \tparam LocMatrix the location matrix type
 */
template <class LocMatrix>
class Newton
{
public:
	/** \brief type of a matrix element */
	typedef double Scalar;

	/** Constructor storing the source and receiver locations */
	Newton(LocMatrix const &xr, LocMatrix const &xs) : xr(xr), xs(xs)
	{
	}

	/** \brief function evaluating the Newton kernel matrix elements */
	double operator()(int i, int j) const
	{
		return 1./(xr.row(i) - xs.row(j)).norm();
	}

private:
	LocMatrix const &xr;	/** \brief receiver locations */
	LocMatrix const &xs;	/** \brief source locations */
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	// read input arguments from Matlab
	NiHu::mex::real_matrix<double> xs(prhs[0]), xr(prhs[1]);
	NiHu::mex::real_matrix<int> source_clusters(prhs[2]), receiver_clusters(prhs[3]), blocks(prhs[4]);
	NiHu::mex::real_matrix<double> input(prhs[5]);

	// allocate Matlab output
	NiHu::mex::real_matrix<double> output(input.rows(), 1, plhs[0]);
	NiHu::mex::real_matrix<int> ranks(blocks.rows(), 1, plhs[1]);
	output.setZero();

	// instantiate Newton kernel matrix
	Newton<NiHu::mex::real_matrix<double> > M(xr, xs);

	ACA aca;
	aca.multiply(M, receiver_clusters, source_clusters, blocks, input, output, 1e-3, 50, ranks);
}

