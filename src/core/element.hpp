// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file element.hpp
 * \ingroup funcspace
 * \brief Declaration of element classes and specialisations
 */
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <iostream>
#include "../util/crtp_base.hpp"
#include "../tmp/bool.hpp"
#include "shapeset.hpp"

/** \brief class describing the overlapping state of two elements */
class element_overlapping
{
public:
	/** \brief return number of coincident nodes */
	unsigned get_num(void) const
	{
		return m_num;
	}

	/** \brief return index of first coincident node in element 1 */
	unsigned get_ind1(void) const
	{
		return m_ind1;
	}

	/** \brief return index of first coincident node in element 2 */
	unsigned get_ind2(void) const
	{
		return m_ind2;
	}

	/** \brief constructor */
	element_overlapping(
		unsigned num = 0, unsigned ind1 = 0, unsigned ind2 = 0) :
		m_num(num), m_ind1(ind1), m_ind2(ind2)
	{
	}

private:
	/** \brief number of coincident nodes */
	unsigned m_num;
	/** \brief start node indices */
	unsigned m_ind1, m_ind2;
};


/** \brief compute surface normal from derivatives in N dimensions
 * \tparam dx_t the derivatives type
 * \tparam N the number of dimensions
 */
template<class scalar, unsigned N>
class normal_impl;

/** \brief specialisation of ::normal_impl for 3D */
template <class scalar>
class normal_impl<scalar, 3>
{
public:
	/**
	 * \brief return the normal vector
	 * \param m the 3x2 derivative matrix
	 * \return the normal vector
	 */
	static Eigen::Matrix<scalar, 3, 1>
	eval(Eigen::Matrix<scalar, 3, 2> const &m)
	{
		return m.col(shape_derivative_index::dXI).cross(m.col(shape_derivative_index::dETA));
	}
};

/** \brief specialisation of ::normal_impl for 2D */
template <class scalar>
class normal_impl<scalar, 2>
{
public:
	/**
	 * \brief return the normal vector
	 * \param m the 2x1 derivative vector
	 * \return the normal vector
	 */
	static Eigen::Matrix<scalar, 2, 1>
	eval(Eigen::Matrix<scalar, 2, 1> const &m)
	{
		return Eigen::Rotation2D<scalar>(-M_PI/2.0) * m;
	}
};


/** \brief compute location derivatives from nodal coordinates */
template <class Derived, unsigned Order>
class location_impl;

/** \brief Traits describing element properties */
namespace element_traits
{
	/** \brief The physical coordinate space of the element */
	template <class Derived>
	struct space_type;

	/** \brief The geometrical shape set of the element */
	template <class Derived>
	struct lset;

	/** \brief Assigns an id to the element type */
	template <class Derived>
	struct id
	{
		enum {
			value = element_traits::space_type<Derived>::type::dimension * 10000 +
			shape_set_traits::id<typename element_traits::lset<Derived>::type>::value
		};
	};

	/** \brief Defines the complexity to determine if the location derivative can be precomputed or not */
	template <class Derived, unsigned Order>
	struct location_complexity : shape_set_traits::shape_complexity<
		typename lset<Derived>::type,
		Order
	> {};

    /** \brief Determines if the normal vector is stored or not */
	template <class Derived>
	struct is_normal_stored : std::integral_constant<bool,
        !std::is_same<
            typename location_complexity<Derived, 1>::type,
            matrix_function_complexity::general
        >::value
    > {};

	/** \brief Matrix that stores the element's corner coordinates */
	template <class Derived>
	struct coords_type
	{
        typedef Eigen::Matrix<
            typename space_type<Derived>::type::scalar_t,
            space_type<Derived>::type::dimension,
            lset<Derived>::type::num_nodes
        > type;
	};

	/** \brief Matrix that stores the physical location's derivatives */
	template <class Derived, unsigned Order>
	struct location_value_type
	{
        typedef Eigen::Matrix<
            typename space_type<Derived>::type::scalar_t,
            space_type<Derived>::type::dimension,
            num_derivatives<
                Order,
                shape_set_traits::domain<
                    typename lset<Derived>::type
                >::type::dimension
            >::value
        > type;
	};

	/** \brief Class that computes or stores the locations */
	template <class Derived, unsigned Order>
	struct location_factory_functor
	{
		typedef conditional_precompute_instance<
			typename location_complexity<Derived, Order>::type,
			location_impl<Derived, Order>,
			typename coords_type<Derived>::type,
			typename shape_set_traits::domain<
                typename lset<Derived>::type
            >::type::xi_t
		> type;
	};

	/** \brief The return type of the physical location's derivatives */
	template <class Derived, unsigned Order>
	struct location_return_type
	{
		typedef typename location_factory_functor<
            Derived,
            Order
        >::type::return_type type;
	};

	template <class Derived>
	struct normal_factory_functor
	{
		typedef conditional_precompute_instance<
			typename location_complexity<Derived, 1>::type,
			normal_impl<
				typename space_type<Derived>::type::scalar_t,
				space_type<Derived>::type::dimension
			>,
			typename location_value_type<Derived, 1>::type
		> type;
	};

	/** \brief The return type of the normal vector */
	template <class Derived>
	struct normal_return_type
	{
		typedef typename normal_factory_functor<
			Derived
		>::type::return_type type;
	};
}

/** \brief metafunction assigning an element to an element tag */
template <class tag>
struct tag2element;

/** \brief type that stores the element's id */
typedef Eigen::Matrix<unsigned, 1, 1> elem_id_t;


template <class Derived, unsigned Order>
class location_impl
{
public:
    typedef typename element_traits::lset<Derived>::type lset_t;
    typedef typename element_traits::coords_type<Derived>::type coords_t;
    typedef typename element_traits::location_value_type<Derived, Order>::type ret_t;
    typedef typename shape_set_traits::domain<lset_t>::type::xi_t xi_t;

    static ret_t eval(coords_t const &coords, xi_t const &xi)
    {
        return coords * lset_t::template eval_shape<Order>(xi);
    }
};


/**
* \brief The geometrical element representation
* \tparam Derived CRTP derived class
*/
template <class Derived>
class element_base
{
public:
	/** \brief self-returning metafunction */
	typedef Derived type;

	/** \brief the element space type */
	typedef typename element_traits::space_type<Derived>::type space_t;
	/** \brief the elements's L-set */
	typedef typename element_traits::lset<Derived>::type lset_t;

	/** \brief the elements's reference domain */
	typedef typename lset_t::domain_t domain_t;
	/** \brief the base domain's scalar variable */
	typedef typename space_t::scalar_t scalar_t;

	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename lset_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef typename lset_t::shape_t L_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef typename lset_t::dshape_t dL_t;

	enum {
		xi_dim = domain_t::dimension,
		x_dim = space_t::dimension,
		/** \brief the element id */
		id = element_traits::id<Derived>::value,
		num_nodes = lset_t::num_nodes,
		num_dd = num_derivatives<2, xi_dim>::value
	};

	/** \brief type of the element's physical location variable \f$ x \f$ */
	typedef typename element_traits::location_value_type<Derived, 0>::type x_t;
	typedef typename element_traits::location_return_type<Derived, 0>::type x_return_type;
	/** \brief type of the gradient of the element's physical location variable \f$ x'_{\xi} \f$ */
	typedef typename element_traits::location_value_type<Derived, 1>::type dx_t;
	typedef typename element_traits::location_return_type<Derived, 1>::type dx_return_type;
	typedef typename element_traits::normal_return_type<Derived>::type normal_return_type;
	/** \brief type of the second derivative of the element's physical location variable \f$ x''_{\xi} \f$ */
	typedef typename element_traits::location_value_type<Derived, 2>::type ddx_t;
	typedef typename element_traits::location_return_type<Derived, 2>::type ddx_return_type;

	/** \brief matrix type that stores the element's corner nodes \f$x_i\f$ */
	typedef Eigen::Matrix<unsigned, num_nodes, 1> nodes_t;
	typedef typename element_traits::coords_type<Derived>::type coords_t;
	/** \brief type that stores the element's id */
	typedef elem_id_t id_t;

protected:
	/** \brief the element's identifier */
	id_t m_id;
	/** \brief the element's nodal indices in the mesh */
	nodes_t m_nodes;
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_t m_coords;
	/** \brief the element's center */
	x_t m_center;
	/** \brief estimated linear size of an element */
	scalar_t m_linear_size_estimate;

	typename element_traits::location_factory_functor<Derived, 0>::type x_computer;
	typename element_traits::location_factory_functor<Derived, 1>::type dx_computer;
	typename element_traits::location_factory_functor<Derived, 2>::type ddx_computer;
	typename element_traits::normal_factory_functor<Derived>::type normal_computer;

	/** \brief the normal return type based on acceleration */
	typedef typename std::conditional<
		element_traits::is_normal_stored<Derived>::value,
		x_t,
		x_t const &
	>::type normal_ret_t;

public:
	NIHU_CRTP_HELPERS

	/** \brief default constructor for std::vector container */
	element_base() {}

	/**
	 * \brief constructor
	 * \param [in] id the elem id
	 * \param [in] nodes the elem node indices
	 * \param [in] coords the elem coordinates
	 */
	element_base(coords_t const &coords, unsigned id = 0, nodes_t const &nodes = nodes_t())
		:m_id((id_t()<<id).finished()),
		m_nodes(nodes),
		m_coords(coords),
		m_center(get_x(domain_t::get_center())),
		x_computer(m_coords, xi_t()),
		dx_computer(m_coords, xi_t()),
		ddx_computer(m_coords, xi_t()),
		normal_computer(get_dx(xi_t()))
	{
	}

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
	 * \param [in] \xi location \f$\xi\f$ in the base domain
	 * \return location \f$x\f$ in the element
	 */
	x_return_type get_x(xi_t const &xi) const
	{
		return x_computer(m_coords, xi);
	}

	/**
	 * \brief return element location gradient
	 * \param [in] xi location \f$\xi\f$ in the base domain
	 * \return location gradient \f$x'_{\xi}\f$ in the element
	 */
	dx_return_type get_dx(xi_t const &xi) const
	{
		return dx_computer(m_coords, xi);
	}

	/**
	 * \brief return element location second derivative matrix
	 * \param [in] xi location \f$ \xi \f$ in the base domain
	 * \return location gradient \f$ x''_{\xi} \f$ in the element
	 */
	ddx_return_type get_ddx(xi_t const &xi) const
	{
		return ddx_computer(m_coords, xi);
	}

	/**
	 * \brief return element center
	 * \return location \f$x_0\f$ in the element center
	 */
	x_t const &get_center(void) const
	{
		return m_center;
	}

	/**
	 * \brief return element normal
	 * \param [in] xi location in the reference domain
	 * \return the element normal
	 */
	normal_return_type get_normal(xi_t const &xi) const
	{
		return normal_computer(dx_computer(m_coords, xi));
	}

	/**
	 * \brief return linear size estimate
	 * \return linear size estimate
	 */
	scalar_t const &get_linear_size_estimate(void) const
	{
		return m_linear_size_estimate;
	}

	/** \brief print debug information on an element to standard output */
	void print_debug(void) const
	{
		std::cout << "Element type id: " << id << std::endl;
		std::cout << "Element id: " << m_id << std::endl;
	}

	/**
	 * \brief check overlapping state with other element
	 * \tparam OtherElem type of the other element
	 * \param [in] other reference to the other element
	 * \return structure containing overlapping information
	 */
	template <class OtherElem>
	element_overlapping get_overlapping(OtherElem const &other) const
	{
		unsigned num_coinc = 0;	// number of coincident nodes
		unsigned start_ind1 = 0, start_ind2 = 0;

		auto const &otherNodes = other.get_nodes();
		for (unsigned i = 0; i < num_nodes; ++i)
		{
			for (unsigned j = 0; j < OtherElem::num_nodes; ++j)
			{
				if (m_nodes[i] == otherNodes[j])
				{
					if (num_coinc == 0)
					{
						start_ind1 = i;
						start_ind2 = j;
					}
					else
					{
						if (i < start_ind1
							|| (i==num_nodes-1
							&& start_ind1 == 0
							&& num_coinc == 1))
							start_ind1 = i;
						if (j < start_ind2 ||
							(j==OtherElem::num_nodes-1 &&
							start_ind2 == 0 &&
							num_coinc == 1))
							start_ind2 = j;
					}
					num_coinc++;
				}
			}
		}
		if (num_coinc == 0)
			return element_overlapping();
		return element_overlapping(num_coinc, start_ind1, start_ind2);
	}
};

/** \brief compare two element id's by value
 * \param [in] lhs the left hand side of comparison
 * \param [in] rhs the right hand side of comparison
 * \return true if lhs < rhs
 */
bool operator<(elem_id_t const &lhs, elem_id_t const &rhs)
{
	return lhs(0,0) < rhs(0,0);
}

// forward declaration
template <class LSet, class scalar_t>
class general_surface_element;

/** \brief specialisation of element_traits for the general surface element */
namespace element_traits
{
	template <class LSet, class Scalar>
	struct space_type<general_surface_element<LSet, Scalar> >
        : space<Scalar, LSet::domain_t::dimension + 1> {};

	template <class LSet, class Scalar>
	struct lset<general_surface_element<LSet, Scalar> >
	{
		typedef LSet type;
	};
}

/** \brief class describing a general surface element that computes its normal in the general way
 * \tparam LSet type of the geometry shape set
 */
template <class LSet, class scalar_t>
class general_surface_element :
	public element_base<general_surface_element<LSet, scalar_t> >
{
public:
	/** \brief the base class type */
	typedef element_base<general_surface_element<LSet, scalar_t> > base_t;

	/** \brief type of the coordinate matrix */
	typedef typename base_t::coords_t coords_t;
	/** \brief type of the node index vector */
	typedef typename base_t::nodes_t nodes_t;
	/** \brief type of the coordinate  on the standard element */
	typedef typename base_t::xi_t xi_t;
	/** \brief type of the coordinate vector */
	typedef typename base_t::x_t x_t;
	/** \brief type of the coordinate derivative vector */
	typedef typename base_t::dx_t dx_t;
	/** \brief type of the reference domain */
	typedef typename base_t::domain_t domain_t;

	/** \brief constructor
	 * \param [in] coords the coordinate matrix
	 * \param [in] id element id
	 * \param [in] nodes the nodal index vector
	 */
	general_surface_element(
		coords_t const &coords,
		unsigned id = 0,
		nodes_t const &nodes = nodes_t())
		: base_t(coords, id, nodes)
	{
		base_t::m_linear_size_estimate = sqrt(base_t::get_normal(domain_t::get_center()).norm() * domain_t::get_volume());
	}
};

#include "../library/line_2_shape_set.hpp"
#include "../library/tria_2_shape_set.hpp"
#include "../library/quad_2_shape_set.hpp"

typedef general_surface_element<line_1_shape_set, space_2d::scalar_t> line_1_elem;
typedef general_surface_element<line_2_shape_set, space_2d::scalar_t> line_2_elem;
typedef general_surface_element<tria_1_shape_set, space_3d::scalar_t> tria_1_elem;
typedef general_surface_element<tria_2_shape_set, space_3d::scalar_t> tria_2_elem;
typedef general_surface_element<quad_1_shape_set, space_3d::scalar_t> quad_1_elem;
typedef general_surface_element<quad_2_shape_set, space_3d::scalar_t> quad_2_elem;


struct line_1_tag {};
template <> struct tag2element<line_1_tag> : line_1_elem {};
struct line_2_tag {};
template <> struct tag2element<line_2_tag> : line_2_elem {};


struct tria_1_tag {};
template <> struct tag2element<tria_1_tag> : tria_1_elem {};
struct tria_2_tag {};
template <> struct tag2element<tria_2_tag> : tria_2_elem {};


struct quad_1_tag {};
template <> struct tag2element<quad_1_tag> : quad_1_elem {};
struct quad_2_tag {};
template <> struct tag2element<quad_2_tag> : quad_2_elem {};

#endif

