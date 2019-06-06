/// \file matrix_free.hpp
/// \brief an Eigen::Matrix adaptor for the fmm::fmm_matrix class 

#ifndef FMM_MATRIX_FREE_HPP_INCLUDED
#define FMM_MATRIX_FREE_HPP_INCLUDED

#include "fmm_matrix.hpp"
#include "../util/real_part_type.hpp"

#include <Eigen/SparseCore>

#include <complex>
#include <type_traits>


namespace NiHu
{
namespace fmm
{
/// \brief an Eigen::Matrix adaptor for the fmm_matrix class
/// \tparam FmmMatrix the fmm_matrix type
template <class FmmMatrix>
class matrix_free;
}
}

namespace Eigen {
namespace internal {
// matrix_free looks-like a SparseMatrix, so let's inherits its traits:
template<class FmmMatrix>
struct traits<NiHu::fmm::matrix_free<FmmMatrix> >
	: public Eigen::internal::traits<Eigen::SparseMatrix<typename FmmMatrix::scalar_t> >
{};
}
}

namespace NiHu
{
namespace fmm
{

/// \brief an Eigen::Matrix adaptor for the fmm_matrix class
/// \tparam FmmMatrix the fmm_matrix type
template <class FmmMatrix>
class matrix_free
	: public Eigen::EigenBase<matrix_free<FmmMatrix> >
{
public:
	typedef FmmMatrix fmm_matrix_t;
	typedef typename fmm_matrix_t::scalar_t scalar_t;

	typedef scalar_t Scalar;
	typedef typename NiHu::real_part_type<Scalar>::type RealScalar;
	typedef int StorageIndex;
	enum {
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> vector_t;

	matrix_free(FmmMatrix &mat)
		: m_mat(mat)
	{
	}

	Eigen::Index rows() const
	{
		return m_mat.rows();
	}

	Eigen::Index cols() const
	{
		return m_mat.cols();
	}

	template <class X>
	typename FmmMatrix::response_t
		operator*(Eigen::MatrixBase<X> const &x) const
	{
		return m_mat * x - x / 2.0;
	}

	vector_t get_diagonal() const
	{
		return m_mat.get_diagonal();
	}

private:
	fmm_matrix_t &m_mat;
};


template <class FmmMatrix>
matrix_free<FmmMatrix>
create_matrix_free(FmmMatrix &fmm_mtx)
{
	return matrix_free<FmmMatrix>(fmm_mtx);
}

} // end of namespace fmm
} // namespace NiHu
#endif // FMM_MATRIX_FREE_HPP_INCLUDED
