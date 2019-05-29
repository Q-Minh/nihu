#ifndef PRECONDITIONER_HPP_INCLUDED
#define PRECONDITIONER_HPP_INCLUDED

namespace NiHu
{
namespace fmm
{

template <class Scalar>
class diagonal_preconditioner
{
public:
	typedef Scalar scalar_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> vector_t;

	diagonal_preconditioner()
	{
	}

	template<typename MatrixType>
	explicit diagonal_preconditioner(MatrixType const &mat)
		: m_inv_diag(mat.cols())
	{
		compute(mat);
	}

	template<typename MatrixType>
	diagonal_preconditioner& analyzePattern(const MatrixType&)
	{
		return *this;
	}

	template<typename MatrixType>
	diagonal_preconditioner& factorize(MatrixType const &mat)
	{
		m_inv_diag = mat.get_diagonal();
		for (size_t i = 0; i < m_inv_diag.rows(); ++i)
			m_inv_diag(i, 0) = 1.0 / m_inv_diag(i, 0);
		return *this;
	}

	template<typename MatrixType>
	diagonal_preconditioner& compute(MatrixType const &mat)
	{
		return factorize(mat);
	}

	template<typename Rhs>
	inline vector_t solve(Rhs const &b) const
	{
		return m_inv_diag.array() * b.array() ;
	}

	Eigen::ComputationInfo info()
	{
		return Eigen::Success;
	}

private:
	vector_t m_inv_diag;
};

} // end of namespace fmm
} // namespace NiHu

#endif
