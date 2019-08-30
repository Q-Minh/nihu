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
	typedef Eigen::Matrix<
		scalar_t,
		Eigen::Dynamic,
		dimension,
		Dim == 1 ? Eigen::ColMajor : Eigen::RowMajor // column vector can not be RowMajor
	> rhs_matrix_t;

	kron_identity()
	{
	}

	explicit kron_identity(Lhs&& lhs)
		: m_lhs(std::forward<Lhs>(lhs))
	{
	}

	Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>
		operator*(Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> const& rhs) const
	{
		/// \todo reimplement this functionality using Eigen::Reshape when Eigen 3.3.9 becomes stable

		// allocate vector result
		Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> res(m_lhs.rows() * dimension, 1);
		// compute result
		Eigen::Map<rhs_matrix_t>(res.data(), m_lhs.rows(), dimension) =
			m_lhs *
			Eigen::Map<const rhs_matrix_t>(rhs.data(), m_lhs.cols(), dimension);
		return res;
	}

private:
	Lhs m_lhs;
};

template <size_t Dim, class Lhs>
kron_identity<Lhs, Dim> create_kron_identity(Lhs&& lhs)
{
	return kron_identity<Lhs, Dim>(std::forward<Lhs>(lhs));
}

} // end of namespace fmm
} // end of namespace NiHu

#endif /* KRON_IDENTITY_HPP_INCLUDED */
