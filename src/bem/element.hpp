/**
* \file element.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of class Element and its specialisations
*/
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include "shapeset.hpp"

template <unsigned nDim>
struct normal_vector
{
	static_assert(nDim <= 2, "Element normal is not yet implemented for this dimension");
};

template <>
struct normal_vector<1>
{
	Eigen::Matrix<double, 1, 2> operator() (Eigen::Matrix<double, 1, 2> const &dx) const
	{
		Eigen::Matrix<double, 1, 2> res;
		res << dx(1), -dx(0);
		return res;
	}
};

template <>
struct normal_vector<2>
{
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
class Element
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
	typedef Matrix<double, 1, x_dim> x_t;
	/** \brief type of the gradient of the element's independent location variable \f$x'_{\xi}\f$ */
	typedef Matrix<double, xi_dim, x_dim> dx_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Matrix<double, num_nodes, x_dim> coords_t;

protected:
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_t coords;

public:
	/**
	* \brief default constructor
	*/
	Element() {}

	/**
	* \brief constructor
	* \param coords location of corners \f$x_i\f$
	*/
	Element(coords_t const &coords) : coords(coords) {}

	/**
	* \brief return element location
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location \f$x\f$ in the element
	*/
	x_t get_x(xi_t const &xi) const
	{
		return lset_t::eval_L(xi).transpose() * coords;
	}

	/**
	* \brief return element location gradient
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location gradient \f$x'_{\xi}\f$ in the element
	*/
	dx_t get_dx(xi_t const &xi) const
	{
		return lset_t::eval_dL(xi).transpose() * coords;
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
typedef Element<line_1_shape_set, 2> line_1_elem;

/** \brief a linear triangle element in 3D space */
typedef Element<tria_1_shape_set, 3> tria_1_elem;
/** \brief a linear parallelogram element in 3D space */
typedef Element<parallelogram_shape_set, 3> parallelogram_elem;
/** \brief a linear quad element in 3D space */
typedef Element<quad_1_shape_set, 3> quad_1_elem;

#endif

