/**
* \file element.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of class Element and its specialisations
*/
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include "shapeset.hpp"

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
	static int const x_dim = Dimension;
	/** \brief the elements's L-set */
	typedef LSet lset;
	/** \brief the dimension of the element's domain variable \f$xi\f$ */
	static int const xi_dim = lset::domain::dimension;

	/** \brief number of shape functions in the set, inherited from the LSet */
	static int const num_nodes = lset::num_nodes;
	/** \brief type of the shape functions' independent variable \f$\xi\f$, inherited from the LSet  */
	typedef typename lset::xi_type xi_type;
	/** \brief type of an \f$L(\xi)\f$ vector, inherited from the LSet */
	typedef typename lset::L_type L_type;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix, inherited from the LSet */
	typedef typename lset::dL_type dL_type;

	/** \brief type of the element's independent location variable \f$x\f$ */
	typedef Matrix<double, 1, x_dim> x_type;
	/** \brief type of the gradient of the element's independent location variable \f$x'_{\xi}\f$ */
	typedef Matrix<double, xi_dim, x_dim> dx_type;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Matrix<double, num_nodes, x_dim> coords_type;

protected:
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_type coords;

public:
	/**
	* \brief constructor
	* \param coords location of corners \f$x_i\f$
	*/
	Element(coords_type const &coords) : coords(coords) {}

	/**
	* \brief return element location
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location \f$x\f$ in the element
	*/
	x_type get_x(xi_type const &xi) const
	{
		return lset::eval_L(xi).transpose() * coords;
	}

	/**
	* \brief return element location gradient
	* \param \xi location \f$\xi\f$ in the base domain
	* \return location gradient \f$x'_{\xi}\f$ in the element
	*/
	dx_type get_dx(xi_type const &xi) const
	{
		return lset::eval_dL(xi).transpose() * coords;
	}

	/**
	* \brief return element normal
	* \param \xi location \f$\xi\f$ in the base domain
	* \return element normal vector \f$n\f$ in the element
	*/
	x_type get_normal(xi_type const &xi) const
	{
		static_assert(xi_dim == x_dim-1, "Element does not have normal");
		static_assert(xi_dim <= 2, "Element normal is not yet implemented for this dimension");
		if (xi_dim == 1) // compile-time switch
		{
				dx_type dx = get_dx(xi);
				dx << dx(1), -dx(0);
				return dx;
		}
		else if (xi_dim == 2)
		{
//			dx_type dx = get_dx(xi);
//			return dx.row(0).cross(dx.row(1));
		}
		else
		{
		}
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

