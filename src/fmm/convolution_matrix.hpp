/**
 * \file convolution_matrix.hpp
 * \brief implemenetation of class template convolution_matrix
 */

#ifndef CONVOLUTION_MATRIX_HPP_INCLUDED
#define CONVOLUTION_MATRIX_HPP_INCLUDED

#include <complex>
#include <Eigen/Dense>

namespace NiHu
{
namespace fmm
{

/** \brief class performing convolution
 * \tparam Scalar the scalar type of the convolution
 */
template <class Scalar>
class convolution_matrix
{
public:
	/** \brief template parameter as nested type */
	typedef Scalar scalar_t;
	/** \brief the column vector type */
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> vector_t;

	/** \brief constructor */
	convolution_matrix()
	{
	}

	/** \brief constructor
	 * \param [in] N the output size of convolution
	 * \param [in] M the input size of convolution
	 * \param [in] diag_coeffs the coefficients of the weight function
	 */
	convolution_matrix(size_t N, size_t M, vector_t const &diag_coeffs)
		: m_N(N)
		, m_M(M)
		, m_diag_coeffs(diag_coeffs)
	{
	}

	/** \brief multiply the convolution matrix by a vector from the right
	 * \param [in] rhs the input of convolution
	 * \return the output of convolution
	 */
	vector_t operator *(vector_t const &rhs) const
	{
		if (rhs.rows() != Eigen::Index(2 * m_M + 1))
			throw std::runtime_error("Invalid input size in convolution matrix");

		vector_t res = vector_t::Zero(2 * m_N + 1, 1);
		int L = int((m_diag_coeffs.rows() - 1) / 2);
		// perform convolution by definition
		for (int n = -int(m_N); n <= m_N; ++n)
		{
			for (int m = -int(m_M); m <= +m_M; ++m)
			{
				int c = n - m;
				if (c >= L || c <= -L)
					continue;
				res(n + m_N) += m_diag_coeffs(c + L) * rhs(m + m_M);
			}
		}
		return res;
	}

private:
	size_t m_N;
	size_t m_M;
	vector_t m_diag_coeffs;
};

} // end of namespace fmm
} // namespace NiHu

#endif // CONVOLUTION_MATRIX_HPP_INCLUDED
