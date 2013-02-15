/**
* \file element.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of class Element and its specialisations
*/
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <Eigen/Dense>
#include <type_traits>

#include "shapeset.hpp"

/**
* \brief functor computing normal vector
* \tparam nDim number of dimensions
*/
template <unsigned nDim>
struct normal_vector
{
	static_assert(nDim <= 2, "Element normal is not yet implemented for this dimension");
};

/**
* \brief specialisation of normal vector for the 1D case
*/
template <>
struct normal_vector<1>
{
	/**
	* \brief compute normal vector in 1D
	* \param dx position derivative matrix
	*/
	Eigen::Matrix<double, 1, 2> operator() (Eigen::Matrix<double, 1, 2> const &dx) const
	{
		Eigen::Matrix<double, 1, 2> res;
		res << dx(1), -dx(0);
		return res;
	}
};

/**
* \brief specialisation of normal vector for the 2D case
*/
template <>
struct normal_vector<2>
{
	/**
	* \brief compute normal vector in two dimensions
	* \param dx position derivative matrix
	*/
	Eigen::Matrix<double, 1, 3> operator() (Eigen::Matrix<double, 2, 3> const &dx) const
	{
		return dx.row(0).cross(dx.row(1));
	}
};


/**
* \brief The geometrical element representation
* \tparam LSet the shape function set describing the geometrical behaviour
* \tparam Dimension the dimensionality of the elements space
* \details The element is defined by its L-set, dimension and the nodal coordinates.
* The class provides a method to compute the location \f$x\f$
*/
template <class LSet, unsigned Dimension>
class element
{
public:
	/** \brief the dimension of the element's location variable \f$x\f$ */
	static unsigned const x_dim = Dimension;
	/** \brief the elements's L-set */
	typedef LSet lset_t;
	/** \brief the elements's domain */
	typedef typename lset_t::domain_t domain_t;
	/** \brief the dimension of the element's domain variable \f$xi\f$ */
	static unsigned const xi_dim = domain_t::dimension;

	/** \brief number of shape functions in the set, inherited from the LSet */
	static unsigned const num_nodes = lset_t::num_nodes;
	/** \brief type of the shape functions' independent variable \f$\xi\f$, inherited from the LSet  */
	typedef typename lset_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector, inherited from the LSet */
	typedef typename lset_t::L_t L_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix, inherited from the LSet */
	typedef typename lset_t::dL_t dL_t;

	/** \brief type of the element's independent location variable \f$x\f$ */
	typedef Eigen::Matrix<double, 1, x_dim> x_t;
	/** \brief type of the gradient of the element's independent location variable \f$x'_{\xi}\f$ */
	typedef Eigen::Matrix<double, xi_dim, x_dim> dx_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Eigen::Matrix<unsigned, 1, num_nodes> nodes_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Eigen::Matrix<double, num_nodes, x_dim> coords_t;
	/** \brief type that stores the element's id */
	typedef Eigen::Matrix<unsigned, 1, 1> id_t;

protected:
	/** \brief the element's identifier */
	id_t m_id;
	/** \brief the element's nodal indices in the mesh */
	nodes_t m_nodes;
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_t m_coords;

public:
	/**
	* \brief default constructor for std::vector container
	*/
	element() {}

	/**
	* \brief constructor
	* \param id the elem id
	* \param nodes the elem node indices
	* \param coords the elem coordinates
	*/
	element(unsigned id, nodes_t const &nodes, coords_t const &coords) : m_nodes(nodes), m_coords(coords)
	{
		this->m_id << id;
	}

	/**
	* \brief constructor
	* \param coords the elem coordinates
	*/
	element(coords_t const &coords) : m_id(id_t()), m_nodes(nodes_t()), m_coords(coords) {}

	/**
	* \brief return elem id
	* \return elem id
	*/
	id_t const &get_id(void) const
	{
		return m_id;
	}

	/**
	* \brief return nodes
	* \return elem nodes
	*/
	nodes_t const &get_nodes(void) const
	{
		return m_nodes;
	}

	/**
	* \brief return coordinates
	* \return elem coordinates
	*/
	coords_t const &get_coords(void) const
	{
		return m_coords;
	}

	/**
	* \brief return element location
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location \f$x\f$ in the element
	*/
	x_t get_x(xi_t const &xi) const
	{
		return lset_t::eval_shape(xi).transpose() * m_coords;
	}

	/**
	* \brief return element center
	* \return location \f$x_0\f$ in the element center
	*/
	x_t get_center(void) const
	{
		return get_x(domain_t::get_center());
	}

	/**
	* \brief return element location gradient
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location gradient \f$x'_{\xi}\f$ in the element
	*/
	dx_t get_dx(xi_t const &xi) const
	{
		return lset_t::eval_dshape(xi).transpose() * m_coords;
	}

	/**
	* \brief return element normal
	* \param \xi location \f$\xi\f$ in the base domain
	* \return element normal vector \f$n\f$ in the element
	*/
	x_t get_normal(xi_t const &xi) const
	{
		static_assert(xi_dim == x_dim-1, "Element does not have normal");
		normal_vector<xi_dim> n;
		return n(get_dx(xi));
	}
};

/** \brief a linear triangle element in 2D space */
typedef element<line_1_shape_set, 2> line_1_elem;

/** \brief a linear triangle element in 3D space */
typedef element<tria_1_shape_set, 3> tria_1_elem;
/** \brief a linear parallelogram element in 3D space */
typedef element<parallelogram_shape_set, 3> parallelogram_elem;
/** \brief a linear quad element in 3D space */
typedef element<quad_1_shape_set, 3> quad_1_elem;

#endif

