#ifndef ACA_HPP_INCLUDED
#define ACA_HPP_INCLUDED

#include <Eigen/Dense>

/** \brief Block matrix representation of a matrix functor */
template <class Matrix, class IndexVector>
class Block
{
public:
	/** \brief the scalar type of the matrix */
	typedef typename std::decay<Matrix>::type::Scalar Scalar;

	/** constructor */
	Block(Matrix M, IndexVector I, IndexVector J)
		: M(std::forward<Matrix>(M)), I(std::forward<IndexVector>(I)), J(std::forward<IndexVector>(J))
	{
	}

	/** index operator */
	Scalar operator()(int i, int j) const 
	{
		return M(I[i], J[j]);
	}

private:
	Matrix M;			/** \brief the matrix */
	IndexVector I, J;	/** \brief the row and column indices */
};


template <class Matrix, class IndexVector>
Block<Matrix, IndexVector>
createBlock(Matrix &&M, IndexVector &&I, IndexVector &&J)
{
	return Block<Matrix, IndexVector>(std::forward<Matrix>(M), std::forward<IndexVector>(I), std::forward<IndexVector>(J));
}


/** \brief perform low rank expansion with a given accuracy and maximum rank */
template <class Matrix, class Result>
int low_rank_approx(Matrix const &M, int nRows, int nCols, double eps, int R, Eigen::DenseBase<Result> &U, Eigen::DenseBase<Result> &V)
{
	// result counter
	int r = 0;

	// estimate typical matrix entry magnitude
	double magest = std::abs(M(nRows/2,nCols/2));

	// start row
	int i = 0;

	// cumulative norm estimate
	double S2 = 0.;

	// iteration counter, R iterations should always be enough
	for (int k = 0; k < R; ++k)
	{
		// compute i-th matrix row
		for (int s = 0; s < nCols; ++s)
			V(s,r) = M(i,s);
		// subtract already approximated parts
		if (r > 0)
			V.col(r) -= V.leftCols(r) * U.block(i,0,1,r).transpose();

		// search row coefficient with maximal magnitude
		int j = 0;
		for (int s = 1; s < nCols; ++s)
			if (std::abs(V(s,r)) > std::abs(V(j,r)))
				j = s;

		// if full zero row has bene found we are ready
		if (std::abs(V(j,r)) / magest < 1e-10)
			break;

		// compute j-th column
		for (int s = 0; s < nRows; ++s)
			U(s,r) = M(s,j);
		if (r > 0)
			U.col(r) -= U.leftCols(r) * V.block(j,0,1,r).transpose();

		// store row and col into result
		U.col(r) /= U(i,r);

		// indicate that r rows and columns are ready
		++r;

		// update estimate of cumulated norm
		S2 = S2 + U.col(r-1).squaredNorm() * V.col(r-1).squaredNorm();
		for (int s = 0; s < r-1; ++s)
			S2 += 2.0 * std::real(U.col(s).dot(U.col(r-1)) * V.col(s).dot(V.col(r-1)));

		// check error
		double err = U.col(r-1).norm() * V.col(r-1).norm() / std::sqrt(S2);
		if (err < eps)
			break;

		// go to next row
		++i;
	}
	
	return r;
}

template <class Matrix, class Input, class Output>
void aca_multiply(Matrix const &M,
	int const *rowIndices, int const *rowLimits, int nRowClusters,
	int const *colIndices, int const *colLimits, int nColClusters,
	int nBlocks, int const *rowBlocks, int const *colBlocks,
	Input const &input, Output &output,
	double eps, int maxRank)
{
	typedef typename Matrix::Scalar Scalar;
	
	// determine maxial row block size
	int nMaxRows = 0;
	for (unsigned i = 0; i < nRowClusters; ++i)
	{
		int nRows = rowLimits[i+1] - rowLimits[i];
		if (nRows > nMaxRows)
			nMaxRows = nRows;
	}
	
	// determine maximal column block size
	int nMaxCols = 0;
	for (unsigned i = 0; i < nColClusters; ++i)
	{
		int nCols = colLimits[i+1] - colLimits[i];
		if (nCols > nMaxCols)
			nMaxCols = nCols;
	}
	
	// allocate memory for U and V
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> U(nMaxRows, maxRank), V(nMaxCols, maxRank);
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Inter(maxRank,1);
		

	// traverse blocks
	for (int iBlock = 0; iBlock < nBlocks; ++iBlock)
	{
		// get row and column cluster indices
		int rowCluster = rowBlocks[iBlock], colCluster = colBlocks[iBlock];
		
		// get cluster size
		int nRows = rowLimits[rowCluster+1]-rowLimits[rowCluster];
		int nCols = colLimits[colCluster+1]-colLimits[colCluster];
		
		// get cluster indices
		int const *rows = rowIndices + rowLimits[rowCluster];
		int const *cols = colIndices + colLimits[colCluster];
		
		auto u = U.topRows(nRows);
		auto v = V.topRows(nCols);
		
		// compute LRA
		int r = low_rank_approx(createBlock(M, rows, cols), nRows, nCols, eps, maxRank, u, v);
		
		auto vv = V.leftCols(r);
		auto uu = U.leftCols(r);
		
		auto inter = Inter.topRows(r);
		inter.setZero();
		for (int j = 0; j < nCols; ++j)
			inter += vv.row(j).transpose() * input[ cols[j] ];
		for (int i = 0; i < nRows; ++i)
			output[ rows[i] ] += uu.row(i)*inter;
	}
}

#endif // ACA_HPP_INCLUDED

