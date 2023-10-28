/** 
 * \file hat_matrix.h
 * \brief implementation of class fmm::hat_matrix
 * \ingroup fmm_util
 */

#ifndef HAT_MATRIX_H_INCLUDED
#define HAT_MATRIX_H_INCLUDED

#include <Eigen/Core>

namespace NiHu
{
namespace fmm
{

class hat_matrix
	: private Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
{
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> base_t;

	base_t const &base() const
	{
		return *(static_cast<base_t const *>(this));
	}

	base_t &base()
	{
		return *(static_cast<base_t *>(this));
	}

public:
	hat_matrix()
	{
	}

	explicit hat_matrix(size_t N)
		: base_t(base_t::Zero(N+1, 2*N+1))
		, m_N(N)
	{
	}

	double &operator()(Eigen::Index n, Eigen::Index m)
	{
		return base()(n + 1, m_N + 1 + m);
	}

	double operator()(Eigen::Index n, Eigen::Index m) const
	{
		return base()(n + 1, m_N + 1 + m);
	}

private:
	size_t m_N;
};

} // end of namespace fmm
} // end of namespace NiHu

#endif /* HAT_MATRIX_H_INCLUDED */
