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
Block<Matrix, IndexVector> createBlock(Matrix &&M, IndexVector &&I, IndexVector &&J)
{
	return Block<Matrix, IndexVector>(std::forward<Matrix>(M), std::forward<IndexVector>(I), std::forward<IndexVector>(J));
}


template <class Matrix, int R>
class ACA
{
public:
	typedef typename Matrix::Scalar Scalar;
	typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> Row;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Col;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, Eigen::Dynamic, R> Utype, Vtype;

	ACA(Matrix const &M, int nRows, int nCols) : M(M), nRows(nRows), nCols(nCols)
	{
		U.resize(nRows, R);
		V.resize(nCols, R);
	}

	/** \brief perform low rank expansion with a given accuracy */
	void process(double eps)
	{
		U.setZero();
		V.setZero();

		// estimate typical matrix entry magnitude
		double magest = std::abs(M(nRows/2,nCols/2));

		// allocate space for a single row and column
		Row row(1,nCols);
		Col col(nRows,1);

		// result counter
		int r = 0;

		// start row
		int i = 0;

		// cumulative norm estimate
		double S2 = 0.;

		// iteration counter, R should be enough
		for (int k = 0; k < R; ++k)
		{
			// compute i-th matrix row
			for (int s = 0; s < nCols; ++s)
				row(s) = M(i,s);
			if (r > 0)
				row -= U.block(i,0,1,r) * V.leftCols(r).transpose();

			// search max abs row coefficient
			int j = 0;
			for (int s = 1; s < nCols; ++s)
				if (std::abs(row(s)) > std::abs(row(j)))
					j = s;

			if (std::abs(row(j)) / magest < 1e-8)
				break;

			// compute j-th column
			for (int s = 0; s < nRows; ++s)
				col(s) = M(s,j);
			if (r > 0)
				col -= U.leftCols(r) * V.block(j,0,1,r).transpose();

			// store row and col into result
			V.col(r) = row.transpose();
			U.col(r) = col/col(i);
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

		// adjust size
		U = (U.leftCols(r)).eval();
		V = (V.leftCols(r)).eval();
	}

	/** \brief return first term of low rank expansion */
	Utype const &get_U(void) const
	{
		return U;
	}

	/** \brief return second term of low rank expansion */
	Vtype const &get_V(void) const
	{
		return V;
	}

private:
	int nRows, nCols;	/** \brief number of matrix rows and columns */
	Matrix const &M;	/** \brief constant reference to the matrix object */
	Utype U;			/** \brief the first term of the low rank expansion */
	Vtype V;			/** \brief the second term of the low rank expansion */
};

