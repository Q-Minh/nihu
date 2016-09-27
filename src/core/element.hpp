// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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
 * \brief Declaration of element classes
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
	
	std::ostream &print_debug(std::ostream &os) const
	{
		return os << "Element overlapping info: " << std::endl
			<< "Number of coincident nodes: " << m_num << std::endl
			<< "Start index on element 1 : " << m_ind1 << " and on element 2 " << m_ind2;
	}

private:
	/** \brief number of coincident nodes */
	unsigned m_num;
	/** \brief start node indices */
	unsigned m_ind1, m_ind2;
};


/** \brief compute location derivatives from nodal coordinates */
template <class Derived, unsigned Order>
class location_impl;

/** \brief Traits describing element properties */
namespace element_traits
{
	/** \brief The element type's textual id */
	template <class Derived>
	struct name
	{
		static const std::string value;
	};

	/** \brief The physical coordinate space of the element */
	template <class Derived>
	struct space_type;

	/** \brief The geometrical shape set of the element */
	template <class Derived>
	struct lset;

	/** \brief indicates if the element is a surface element or not */
	template <class Derived>
	struct is_surface_element;

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
			num_derivatives(Order, shape_set_traits::domain<typename lset<Derived>::type>::type::dimension)
        > type;
	};

	/** \brief Class that computes or stores the locations */
	template <class Derived, unsigned Order>
	struct location_factory_functor : conditional_precompute_instance<
		typename location_complexity<Derived, Order>::type,
		location_impl<Derived, Order>,
		typename coords_type<Derived>::type,
		typename shape_set_traits::domain<
			typename lset<Derived>::type
		>::type::xi_t
	> {};

	/** \brief The return type of the physical location's derivatives */
	template <class Derived, unsigned Order>
	struct location_return_type
	{
		typedef typename location_factory_functor<
            Derived,
            Order
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
	/** \brief the element's L-set */
	typedef typename element_traits::lset<Derived>::type lset_t;
	/** \brief the element's scalar variable */
	typedef typename space_t::scalar_t scalar_t;

	/** \brief the element's reference domain */
	typedef typename lset_t::domain_t domain_t;
	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename domain_t::xi_t xi_t;

	enum {
		/** \brief the element id */
		id = element_traits::id<Derived>::value,
		/** \brief the number of nodes */
		num_nodes = lset_t::num_nodes,
	};

	/** \brief type of the element's physical location variable \f$ x \f$ */
	typedef typename element_traits::location_value_type<Derived, 0>::type x_t;
	/** \brief return type of the element physical location variable when obtained by get_x */
	typedef typename element_traits::location_return_type<Derived, 0>::type x_return_type;
	/** \brief type of the gradient of the element's physical location variable \f$ x'_{\xi} \f$ */
	typedef typename element_traits::location_value_type<Derived, 1>::type dx_t;
	/** \brief return type of the element physical location derivative */
	typedef typename element_traits::location_return_type<Derived, 1>::type dx_return_type;
	/** \brief type of the second derivative of the element's physical location variable \f$ x''_{\xi} \f$ */
	typedef typename element_traits::location_value_type<Derived, 2>::type ddx_t;
	/** \brief return type of the element physical location second derivative */
	typedef typename element_traits::location_return_type<Derived, 2>::type ddx_return_type;

	/** \brief vector type that stores the element's corner node indices */
	typedef Eigen::Matrix<unsigned, num_nodes, 1> nodes_t;
	/** \brief matrix type that stores the element's corner nodes \f$x_i\f$ */
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

	/** \brief functor computing the locations */
	typename element_traits::location_factory_functor<Derived, 0>::type x_computer;
	/** \brief functor computing the location derivative vector */
	typename element_traits::location_factory_functor<Derived, 1>::type dx_computer;
	/** \brief functor computing the location second derivative matrix */
	typename element_traits::location_factory_functor<Derived, 2>::type ddx_computer;

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
		ddx_computer(m_coords, xi_t())
	{
		// compute linear size estimate (average edge length)
		scalar_t d = 0.0;
		for (unsigned e = 0; e < domain_t::num_edges; ++e)
		{
			auto const &edge = domain_t::get_edge(e);
			x_t x0 = get_x(domain_t::get_corner(edge[0]));
			x_t x1 = get_x(domain_t::get_corner(edge[1]));
			d += (x1-x0).norm();
		}
		d /= domain_t::num_edges;
		m_linear_size_estimate = d;
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
	 * \brief return linear size estimate
	 * \return linear size estimate
	 */
	scalar_t const &get_linear_size_estimate(void) const
	{
		return m_linear_size_estimate;
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
inline bool operator<(elem_id_t const &lhs, elem_id_t const &rhs)
{
	return lhs(0,0) < rhs(0,0);
}


// forward declaration
template <class LSet, class scalar_t>
class surface_element;

// forward declaration
template <class LSet, class scalar_t>
class volume_element;


/** \brief specialisation of element_traits for the volume element */
namespace element_traits
{
	/* the space dimension is computed from the Lset domain's dimension */
	template <class LSet, class Scalar>
	struct space_type<volume_element<LSet, Scalar> >
		: space<Scalar, LSet::domain_t::dimension> {};

	template <class LSet, class Scalar>
	struct lset<volume_element<LSet, Scalar> >
	{
		typedef LSet type;
	};

	template <class LSet, class Scalar>
	struct is_surface_element<volume_element<LSet, Scalar> > : std::false_type {};
}



/** \brief compute surface normal from location derivatives */
template<class Derived, class enable = void>
class normal_impl;

/** \brief specialisation of element_traits for the surface element */
namespace element_traits
{
	/* the space dimension is computed from the Lset domain's dimension (dL+1) */
	template <class LSet, class Scalar>
	struct space_type<surface_element<LSet, Scalar> >
		: space<Scalar, LSet::domain_t::dimension + 1> {};

	template <class LSet, class Scalar>
	struct lset<surface_element<LSet, Scalar> >
	{
		typedef LSet type;
	};

	template <class LSet, class Scalar>
	struct is_surface_element<surface_element<LSet, Scalar> > : std::true_type {};
}

/** \brief specialisation of ::normal_impl for 3D */
template <class Derived>
class normal_impl<Derived, typename std::enable_if<
	element_traits::space_type<Derived>::type::dimension == 3
>::type>
{
	typedef typename element_traits::location_value_type<Derived, 0>::type x_t;
	typedef typename element_traits::location_value_type<Derived, 1>::type dx_t;
public:
	static x_t eval(dx_t const &m)
	{
		return m.col(shape_derivative_index::dXI).cross(m.col(shape_derivative_index::dETA));
	}
};

/** \brief specialisation of ::normal_impl for 2D */
template <class Derived>
class normal_impl<Derived, typename std::enable_if<
	element_traits::space_type<Derived>::type::dimension == 2
>::type>
{
	typedef typename element_traits::location_value_type<Derived, 0>::type x_t;
	typedef typename element_traits::location_value_type<Derived, 1>::type dx_t;
public:
	static x_t eval(dx_t const &m)
	{
		return Eigen::Rotation2D<typename Derived::scalar_t>(-M_PI/2.0) * m;
	}
};

namespace element_traits
{
	/** \brief Class that computes or stores the normals */
	template <class Derived>
	struct normal_factory_functor : conditional_precompute_instance<
		typename location_complexity<Derived, 1>::type,
		normal_impl<Derived>,
		typename location_value_type<Derived, 1>::type
	> {};

	/** \brief The return type of the normal vector */
	template <class Derived>
	struct normal_return_type
	{
		typedef typename normal_factory_functor<
			Derived
		>::type::return_type type;
	};
}

/** \brief class describing a surface element that provides a normal vector
 * \tparam LSet type of the geometry shape set
 */
template <class LSet, class scalar_t>
class surface_element :
	public element_base<surface_element<LSet, scalar_t> >
{
public:
	/** \brief the base class type */
	typedef element_base<surface_element<LSet, scalar_t> > base_t;

	using typename base_t::coords_t;
	using typename base_t::nodes_t;
	using typename base_t::xi_t;
	using typename base_t::x_t;
	using typename base_t::dx_t;

	/** \brief the domain type */
	using domain_t = typename base_t::domain_t;

	using base_t::get_dx;

	/** \brief the return type of function get_normal */
	typedef typename element_traits::normal_return_type<surface_element>::type normal_return_type;

protected:
	/** \brief functor computing the normal vector */
	typename element_traits::normal_factory_functor<surface_element>::type normal_computer;

public:
	/** \brief constructor
	 * \param [in] coords the coordinate matrix
	 * \param [in] id element id
	 * \param [in] nodes the nodal index vector
	 * \todo this sqrt is only valid for 3D surfaces
	 */
	surface_element(
		coords_t const &coords,
		unsigned id = 0,
		nodes_t const &nodes = nodes_t())
		: base_t(coords, id, nodes),
		normal_computer(get_dx(xi_t()))
	{
	}
	
	void flip_inplace(void)
	{
		*this = flip();
	}
	
	surface_element flip(void) const
	{
		return surface_element(this->get_coords().colwise().reverse(), this->get_id()(0), this->get_nodes().reverse());
	}

	/**
	 * \brief return element normal
	 * \param [in] xi location in the reference domain
	 * \return the element normal
	 */
	normal_return_type get_normal(xi_t const &xi = domain_t::get_center()) const
	{
		return normal_computer(get_dx(xi));
	}
};



/** \brief class describing a volume element that has no normal vector
 * \tparam LSet type of the geometry shape set
 */
template <class LSet, class scalar_t>
class volume_element :
	public element_base<volume_element<LSet, scalar_t> >
{
public:
	/** \brief the base class type */
	typedef element_base<volume_element<LSet, scalar_t> > base_t;

	using typename base_t::coords_t;
	using typename base_t::nodes_t;
	using typename base_t::xi_t;
	using typename base_t::x_t;
	using typename base_t::dx_t;

	/** \brief the domain type */
	using domain_t = typename base_t::domain_t;

	using base_t::get_dx;

	/** \brief constructor
	 * \param [in] coords the coordinate matrix
	 * \param [in] id element id
	 * \param [in] nodes the nodal index vector
	 * \todo this linear size estimate is incorrect
	 */
	volume_element(
		coords_t const &coords,
		unsigned id = 0,
		nodes_t const &nodes = nodes_t())
		: base_t(coords, id, nodes)
	{
	}
};


#endif

