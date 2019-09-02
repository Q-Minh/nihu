/** 
 * @file chebyshev_cluster.hpp
 * @brief Implementation of class @ref chebyshev_cluster
 * @ingroup bbfmm
 */

#ifndef CHEBYSHEV_CLUSTER_HPP_INCLUDED
#define CHEBYSHEV_CLUSTER_HPP_INCLUDED

#include "cluster.hpp"
#include "util/misc.hpp"
#include "nd_cheb.hpp"

namespace NiHu
{
namespace fmm
{

/** 
 * @brief Cluster class of the black box fmm
 * @tparam Dim Dimension of the system (location)
 * @tparam Scalar Kernel value's scalar type
 * @tparam FieldDim Dimension of the field variable 
 */
template <size_t Dim, class Scalar, size_t FieldDim>
class chebyshev_cluster;

/** 
 * @brief Traits class of the @ref chebyshev_cluster
 * @tparam Dim Space dimension
 * @tparam Scalar Scalar type
 * @tparam FieldDim Field's dimension
 */
template <size_t Dim, class Scalar, size_t FieldDim>
struct cluster_traits<chebyshev_cluster<Dim, Scalar, FieldDim> >
{
	/// \brief the space's dimension
	static size_t const dimension = Dim;
	/// \brief the multipole contribution's type
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> multipole_t;
	/// \brief the local contribution's type
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> local_t;
};

template <size_t Dim, class Scalar, size_t FieldDim = 1>
class chebyshev_cluster
	: public fmm::cluster_base<chebyshev_cluster<Dim, Scalar> >
{
public:
	/// \brief the CRTP base class
	typedef fmm::cluster_base<chebyshev_cluster<Dim, Scalar> > base_t;
	/// \brief the space's dimension
	static size_t const dimension = base_t::dimension;
	/// \brief the field dimension
	static size_t const field_dimension = FieldDim;
	/// \brief the multipole type
	typedef typename base_t::multipole_t multipole_t;
	/// \brief the local type
	typedef typename base_t::local_t local_t;

	/// \brief type to store Chebyshev nodes
	typedef Eigen::Matrix<double, dimension, Eigen::Dynamic> cheb_nodes_t;

	/// \brief set the Chebyshev order
	/// \param [in] order the Chebyshev order
	void set_chebyshev_order(size_t order)
	{
		m_chebyshev_order = order;
		m_cheb_nodes = chebnodes<dimension>(m_chebyshev_order, this->get_bounding_box());
	}

	/// \brief return the number of elements in the multipole coefficients
	/// \return number of multipole coefficients
	size_t get_data_size() const
	{
		return Ipow(m_chebyshev_order, dimension) * field_dimension;
	}

	/// \brief return zero multipole contribution
	/// \return zero multipole contribution
	multipole_t zero_multipole() const
	{
		return multipole_t::Zero(get_data_size(), 1);
	}

	/// \brief return zero local contribution
	/// \return zero local contribution
	local_t zero_local() const
	{
		return local_t::Zero(get_data_size(), 1);
	}

	/// \brief return the Chebyshev nodes
	/// \return Chebyshev nodes
	cheb_nodes_t const &get_chebyshev_nodes() const
	{
		return m_cheb_nodes;
	}

	/// \brief return the Chebyshev order
	/// \return Chebyshev order
	size_t get_chebyshev_order() const
	{
		return m_chebyshev_order;
	}

private:
	size_t m_chebyshev_order;
	cheb_nodes_t m_cheb_nodes;
};

} // end of namespace fmm
} // end of namespace NiHu

#endif /* CHEBYSHEV_CLUSTER_HPP_INCLUDED */
