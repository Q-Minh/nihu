/** \file aca.hpp
* \details implementation of Adaptive Cross Approximation (ACA)
*/

#ifndef ACA_HPP_INCLUDED
#define ACA_HPP_INCLUDED

#include <Eigen/Dense>
#include <type_traits> // std::decay
#include <utility> // std::forward
#include <vector> // std::forward

template <class Scalar>
class LowRank
{
public:
	LowRank() {}

	template <class UType, class VType>
	LowRank(int row0, int col0, Eigen::MatrixBase<UType> const &U, Eigen::DenseBase<VType> const &V)
		: row0(row0), col0(col0), U(U), V(V) {}

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> const &getU(void) const { return U; }
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> const &getV(void) const { return V; }

	int getRow0(void) const { return row0; } 
	int getCol0(void) const { return col0; } 

private:
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> U;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> V;
	int row0, col0;
};

/** \brief class performing Adaptive Cross Approximation */
class ACA
{
private:
	/** \brief Block matrix representation of a matrix functor */
	template <class Matrix>
	class Block
	{
	public:
		/** \brief the scalar type of the matrix */
		typedef typename std::decay<Matrix>::type::Scalar Scalar;

		/** constructor */
		Block(Matrix M, int row0, int col0)
			: M(std::forward<Matrix>(M)), row0(row0), col0(col0) { }

		/** index operator */
		Scalar operator()(int i, int j) const { return M(row0+i, col0+j); }

	private:
		Matrix M;			/** \brief the matrix */
		int row0, col0;		/** \brief the row and column offsets */
	};


	template <class Matrix>
	static Block<Matrix>
		createBlock(Matrix &&M, int row0, int col0)
	{
		return Block<Matrix>(std::forward<Matrix>(M), row0, col0);
	}


	/** \brief compute low rank approximation of a matrix with a prescribed accuracy and maximum rank
	* \details The low rank approximation is of the form M = U * V.transpose()
	* \tparam Matrix the matrix class
	* \tparam Result the type of the result (Eigen matrix)
	* \param [in] M the matrix to approximate
	* \param [in] nRows the number of matrix rows
	* \param [in] nCols the number of matrix cols
	* \param [in] eps the prescribed approximation accuracy
	* \param [in] R the maximum rank (number of iterations)
	* \param [out] U lhs term of the approximation
	* \param [out] V rhs term of the approximation
	* \return the actual rank of the approximation ( <= R )
	*/
	template <class Matrix, class Result>
	static int low_rank_approx(Matrix const &M, int nRows, int nCols, double eps, int R,
		Eigen::DenseBase<Result> &U, Eigen::DenseBase<Result> &V)
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

			// if full zero row has been found we are ready
			if (std::abs(V(j,r)) / magest < 1e-10)
				break;

			// compute j-th column
			for (int s = 0; s < nRows; ++s)
				U(s,r) = M(s,j);
			// subtract already approximated parts
			if (r > 0)
				U.col(r) -= U.leftCols(r) * V.block(j,0,1,r).transpose();
			// normalise
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

public:
	/** \brief Compute matrix-vector product with multilevel low rank approximation
	* \tparam Matrix	the matrix class
	* \tparam RowArray	type of the array storing the row cluster tree
	* \tparam ColumnArray	type of the array storing the column cluster tree
	* \tparam BlockArray	type of the array storing the block tree
	* \tparam Input	type of the input vector
	* \tparam Output	type of the output vector
	* \param [in] rowClusters	the row cluster tree as Eigen matrix
	* \param [in] colClusters	the column cluster tree as Eigen matrix
	* \param [in] blocks	the block tree as Eigen matrix
	* \param [in] input	the input vector
	* \param [out] output	the output vector
	* \param [in] eps	the ACA approximation error
	* \param [in] maxRank	the maximal ACA approximation rank
	*/
	template <class Matrix, class RowArray, class ColumnArray, class BlockArray, class Input, class Output>
	void multiply(Matrix const &M,
		Eigen::DenseBase<RowArray> const &rowClusters, Eigen::DenseBase<ColumnArray> const &colClusters,
		Eigen::DenseBase<BlockArray> const &blocks,
		Input const &input, Output &output,
		double eps, int maxRank)
	{
		typedef typename Matrix::Scalar Scalar;

		// determine maximal row and column block size for preallocation
		int nMaxRows = rowClusters.col(1).maxCoeff(), nMaxCols = colClusters.col(1).maxCoeff();

		// allocate memory for U, V and internal result vector of LRA
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> U(nMaxRows, maxRank), V(nMaxCols, maxRank);
		Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Z(maxRank,1);

		// track number of matrix element evaluations
		m_nEval = m_matrixSize = 0;

		// traverse blocks
		for (int iBlock = 0; iBlock < blocks.rows(); ++iBlock)
		{
			auto rowCluster(rowClusters.row(blocks(iBlock,0)));	// reference object
			auto colCluster(colClusters.row(blocks(iBlock,1)));	// reference object

			// get cluster size
			int row0 = rowCluster(0), nRows = rowCluster(1);
			int col0 = colCluster(0), nCols = colCluster(1);

			auto u(U.topRows(nRows));	// reference objects
			auto v(V.topRows(nCols));

			// compute low rank decomposition
			int r = low_rank_approx(createBlock(M, row0, col0), nRows, nCols, eps, maxRank, u, v);

			m_nEval += (nRows+nCols)*r;
			m_matrixSize += nRows*nCols;

			// perform multiplication
			auto vv(V.leftCols(r));		// reference objects
			auto uu(U.leftCols(r));
			auto z(Z.topRows(r));

			// z = V' * y
			z.setZero();
			for (int j = 0; j < nCols; ++j)
				z += vv.row(j).transpose() * input( col0+j , 0);
			// x += U * z
			for (int i = 0; i < nRows; ++i)
				output( row0+i, 0 ) += (uu.row(i)*z)(0,0);
		}
	}

	template <class Matrix, class RowArray, class ColumnArray, class BlockArray>
	static std::vector<LowRank<typename Matrix::Scalar> > decompose(Matrix const &M,
		Eigen::DenseBase<RowArray> const &rowClusters, Eigen::DenseBase<ColumnArray> const &colClusters,
		Eigen::DenseBase<BlockArray> const &blocks,
		double eps, int maxRank)
	{
		typedef typename Matrix::Scalar Scalar;
		typedef std::vector<LowRank<Scalar> > ReturnType;

		ReturnType ret;
		ret.reserve(blocks.rows());

		// determine maximal row and column block size for preallocation
		int nMaxRows = rowClusters.col(1).maxCoeff(), nMaxCols = colClusters.col(1).maxCoeff();

		// allocate memory for U, V and internal result vector of LRA
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> U(nMaxRows, maxRank), V(nMaxCols, maxRank);
		Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Z(maxRank,1);

		// traverse blocks
		for (int iBlock = 0; iBlock < blocks.rows(); ++iBlock)
		{
			auto rowCluster(rowClusters.row(blocks(iBlock,0)));	// reference object
			auto colCluster(colClusters.row(blocks(iBlock,1)));	// reference object

			// get cluster size
			int row0 = rowCluster(0), nRows = rowCluster(1);
			int col0 = colCluster(0), nCols = colCluster(1);

			auto u(U.topRows(nRows));	// reference objects
			auto v(V.topRows(nCols));

			// compute low rank decomposition
			int r = low_rank_approx(createBlock(M, row0, col0), nRows, nCols, eps, maxRank, u, v);

			ret.push_back(LowRank<Scalar>(row0, col0, U.leftCols(r), V.leftCols(r)));
		}

		return ret;
	}

public:
	/** \brief return number of matrix evaluations */
	int get_nEval(void) const { return m_nEval; }
	/** \brief number of matrix elements in processed blocks */
	int get_matrixSize(void) const { return m_matrixSize; }
	/** \brief return the size compression ratio */
	double get_sizeCompression(void) const { return double(get_nEval()) / get_matrixSize(); }

private:
	int m_nEval;
	int m_matrixSize;
};

#endif // ACA_HPP_INCLUDED

