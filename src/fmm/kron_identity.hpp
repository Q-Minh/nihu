/**
 * \file kron_identity.hpp
 * \brief kronecker product of a matrix by an Identity matrix
 * \ingroup fmm_util
 */

#ifndef KRON_IDENTITY_HPP_INCLUDED
#define KRON_IDENTITY_HPP_INCLUDED

#include <Eigen/Dense>

#include <type_traits>
#include <utility>

#include "mex.h"

namespace NiHu
{
namespace fmm
{

template <class Lhs, size_t Dim>
class kron_identity
{
public:
	typedef typename std::decay<Lhs>::type lhs_t;
	static unsigned const dimension = Dim;
	typedef typename lhs_t::Scalar scalar_t;
	static Eigen::Index const lhs_rows_compile_time = lhs_t::RowsAtCompileTime;
	static Eigen::Index const result_rows_compile_time =
		lhs_rows_compile_time == Eigen::Dynamic ? Eigen::Dynamic : dimension * lhs_rows_compile_time;

	kron_identity()
	{
	}

	explicit kron_identity(Lhs &&lhs)
		: m_lhs(std::forward<Lhs>(lhs))
	{
	}

	template <class RhsDerived>
	Eigen::Matrix<scalar_t, result_rows_compile_time, RhsDerived::ColsAtCompileTime>
		operator*(Eigen::MatrixBase<RhsDerived> const &rhs) const
	{
		// allocate  result
		Eigen::Matrix<scalar_t, result_rows_compile_time, RhsDerived::ColsAtCompileTime> res(m_lhs.rows() * dimension, rhs.cols());
		res.setZero();
		for (Eigen::Index i = 0; i < m_lhs.rows(); ++i)
			for (Eigen::Index j = 0; j < m_lhs.cols(); ++j)
				res.block(i * dimension, 0, dimension, rhs.cols()) +=
				m_lhs(i, j) * rhs.block(j * dimension, 0, dimension, rhs.cols());
		return res;
	}

private:
	Lhs m_lhs;
};

template <size_t Dim, class Lhs>
kron_identity<Lhs, Dim> create_kron_identity(Lhs &&lhs)
{
	return kron_identity<Lhs, Dim>(std::forward<Lhs>(lhs));
}

} // end of namespace fmm
} // end of namespace NiHu

#endif /* KRON_IDENTITY_HPP_INCLUDED */
