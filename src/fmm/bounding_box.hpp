/** 
 * @file bounding_box.hpp
 * @brief Implementation of class @ref NiHu::fmm::bounding_box
 * @ingroup fmm_clus
 */

#ifndef BOUNDING_BOX_HPP_INCLUDED
#define BOUNDING_BOX_HPP_INCLUDED

#include <Eigen/Dense>

#include <algorithm>
#include <iostream>
#include <cmath>

namespace NiHu
{
namespace fmm
{
/** \brief multidimensional square bounding box
 * \tparam Dim the space dimension
 * \tparam Scalar scalar type of coordinates
 *
 * A bounding_box is a squared box aligned to the coordinate directions,
 * defined by a center and a diameter. bounding_box is the main building block
 * component of the fmm cluster and the cluster_tree.
 */
template <size_t Dim, class Scalar = double>
class bounding_box
{
public:
	/** \brief template parameter as nested constant */
	static size_t const dimension = Dim;
	/** \brief template argument as nested type */
	typedef Scalar scalar_t;
	/** \brief the location type in the bounding box */
	typedef Eigen::Matrix<scalar_t, dimension, 1> location_t;

	/** \brief constructor from center and diameter
	 * \param [in] center the box center
	 * \param [in] diameter the box diameter (edge length)
	 */
	bounding_box(
		location_t const &center = location_t::Zero(),
		scalar_t diameter = 2.0
	)
		: m_center(center)
		, m_diameter(diameter)
	{
	}

	/** \brief constructor from a set of contained nodes
	 * \tparam NodesDerived Eigen derived class of Nodes
	 * \param [in] nodes matrix containing the enclosed nodes in its columns
	 *
	 * The bounding_box is composed so that tightly wraps all the nodes.
	 */
	template <class NodesDerived>
	explicit bounding_box(Eigen::DenseBase<NodesDerived> const &nodes)
		: m_diameter(0.)
	{
		for (size_t d = 0; d < dimension; ++d)
		{
			scalar_t a = nodes.row(d).minCoeff();
			scalar_t b = nodes.row(d).maxCoeff();
			m_center(d) = (a + b) / 2.;
			m_diameter = std::max(m_diameter, b-a);
		}
	}

	/** \brief return center
	 * \return center
	 */
	location_t const &get_center(void) const
	{
		return m_center;
	}

	/** \brief set the center
	 * \param [in] c the new center
	 */
	void set_center(location_t const &c)
	{
		m_center = c;
	}

	/** \brief return diameter
	 * \return diameter
	 */
	scalar_t get_diameter(void) const
	{
		return m_diameter;
	}

	/** \brief set the diameter
	 * \param [in] d the new diameter
	 */
	void set_diameter(scalar_t const &d)
	{
		m_diameter = d;
	}

	/** \brief get a child box
	 * \param [in] idx the child index
	 * \return the specified child box
	 */
	bounding_box get_child(size_t idx) const
	{
		location_t offset = idx2dist(idx) * (m_diameter / 4.0);
		return bounding_box(m_center + offset, m_diameter / 2.0);
	}

	/** \brief determine if a box is adjacent
	 * \param [in] other an other bounding_box
	 * \param [in] tol the relative tolerance of distance measurement
	 * \param [in] (optional) tolerance for distance checking
	 * \return true if the two are adjacent
	 *
	 * Two boxes are adjacent if the distance of their centres is not larger than
	 * the sum of the diameters. This condition is checked for each coordinate
	 * direction.
	 */
	bool is_adjacent(bounding_box const &other, double tol = 1e-3) const
	{
		location_t dist = other.get_center() - this->get_center();
		scalar_t D = (this->get_diameter() + other.get_diameter()) / 2.;
		for (size_t d = 0; d < dimension; ++d)
			if (std::abs(dist(d)) > D * (1.0 + tol))
				return false;
		return true;
	}

	/** \brief convert parent-to-child direction to child index
	 * \param [in] child location of the child
	 * \param [in] parent location of the parent
	 * \return child index
	 */
	static unsigned dist2idx(location_t const &child, location_t const &parent)
	{
		unsigned idx = 0;
		for (size_t d = 0; d < dimension; ++d)
			if (child(d) > parent(d))
				idx |= (1 << d);
		return idx;
	}

	/** \brief convert child index to parent-to-child direction
	 * \param [in] idx the child index
	 * \return parent-to-child direction vector
	 */
	static location_t idx2dist(size_t idx)
	{
		if (idx >= (1 << dimension))
			throw std::invalid_argument("bounding_box child index too large");

		location_t dst = location_t::Zero();
		for (unsigned d = 0; d < dimension; ++d)
			dst(d) = ((idx >> d) & 1) == 1 ? 1.0 : -1.0;
		return dst;
	}

	/** \brief print debug information to an output stream
	 * \param [in] os the output stream
	 */
	void print_debug(std::ostream &os = std::cout) const
	{
		os << "(" << this->get_center().transpose() << ") " << this->get_diameter();
	}

private:
	location_t m_center;
	scalar_t m_diameter;
};

/** \brief inserter operator of a bounding box
 * \tparam Dim the dimension
 * \tparam T the scalar type
 * \param [in, out] os the output stream
 * \param [in] bb the bounding box
 * \return the modified output stream
 */
template <size_t Dim, class T>
std::ostream &operator<<(std::ostream &os, bounding_box<Dim, T> const &bb)
{
	bb.print_debug(os);
	return os;
}
} // namespace fmm
}
#endif // BOUNDING_BOX_HPP_INCLUDED
