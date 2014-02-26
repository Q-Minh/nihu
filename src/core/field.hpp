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
 * \file field.hpp
 * \ingroup funcspace
 * \brief implementation of fields and field views
 * \details Fields are elements extended with a shape set and dof description.
 * Fields can be generated ,,by hand'' or by extending an element with automatic
 * field generation options. This file implements class ::field, class ::field_view
 * and class ::dirac_field. All field-type classes are derived from the CRTP base ::field_base
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "element.hpp"

/** \brief collect options used to convert an element into a field view */
namespace field_option
{
	/** \brief tag to describe an isoparametric field */
	struct isoparametric {};
	/** \brief tag to describe a constant field */
	struct constant {};
}

typedef std::integral_constant<unsigned, 1> _1d;
typedef std::integral_constant<unsigned, 2> _2d;
typedef std::integral_constant<unsigned, 3> _3d;

/** \brief assign a field to a tag
 * \tparam field_tag the tag
 */
template <class field_tag>
struct tag2field;


/** \brief traits class of all fields */
template <class Derived>
struct field_traits;


/** \brief metafunction assigning a default id to a field
* \tparam field_t the field type
*/
template <class field_t>
struct field_id
{
	/** \brief the default field id */
	static unsigned const value =
		element_traits::id<typename field_traits<field_t>::elem_t>::value * 100 +
		field_traits<field_t>::nset_t::num_nodes;
};


/**
 * \brief CRTP base class of all fields
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class field_base
{
public:
	NIHU_CRTP_HELPERS

	/** \brief self returning metafunction */
	typedef Derived type;

	/** \brief the traits class */
	typedef field_traits<Derived> traits_t;
	/** \brief the element type */
	typedef typename traits_t::elem_t elem_t;
	/** \brief the nset type */
	typedef typename traits_t::nset_t nset_t;
	/** \brief the dofs vector type */
	typedef typename traits_t::dofs_t dofs_t;

	enum {
		/** \brief the number of dofs */
		num_dofs = nset_t::num_nodes,
		/** \brief the field id */
		id = field_id<Derived>::value,
		/** \brief the quantity's dimensionality */
		quantity_dimension = traits_t::quantity_dimension
	};

	/** \brief return underlying element */
	elem_t const &get_elem(void) const
	{
		return derived().get_elem();
	}

	/** \brief return DOF vector */
	dofs_t const &get_dofs(void) const
	{
		return derived().get_dofs();
	}
};



/** \brief imlementation class of a general field */
template <class Derived>
class field_impl;

// forward declaration
template <class Field>
class dirac_field;

/** \brief specialisation of ::field_traits for a ::dirac_field
 * \tparam Field the original field that is converted into a dirac_field
 */
template <class Field>
struct field_traits<dirac_field<Field> >
{
	/** \brief the element type */
	typedef typename field_traits<Field>::elem_t elem_t;
	/** \brief the nset type */
	typedef typename field_traits<Field>::nset_t nset_t;
	enum {
		quantity_dimension = field_traits<Field>::quantity_dimension
	};
	/** \brief the dof vector type */
	typedef typename field_traits<Field>::dofs_t dofs_t;
	/** \brief indicates that the field is dirac */
	static bool const is_dirac = true;
};


/** \brief dirac view of a field
 * \tparam Field the original field that is converted into a Dirac field
 */
template <class Field>
class dirac_field :
	public field_base<dirac_field<Field> >,
	public field_impl<Field>
{
public:
	/** \brief self-returning metafunction */
	typedef dirac_field type;

	/** \brief the implementation class type */
	typedef field_impl<Field> impl_t;

	/** \brief shorthand for the element type */
	typedef typename field_traits<dirac_field<Field> >::elem_t elem_t;
	/** \brief store the original field type for the iterator */
	typedef Field field_t;

	using impl_t::get_elem;
	using impl_t::get_dofs;
};


// forward declaration
template <class ElemType, class FieldOption, class Dimension>
class field_view;

/** \brief Specialisation of field_traits for the isoparametric field view case */
template <class ElemType, class Dimension>
struct field_traits<field_view<ElemType, field_option::isoparametric, Dimension> >
{
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief the dof vector type */
	typedef typename elem_t::lset_t nset_t;

	enum { quantity_dimension = Dimension::value };
	/** \brief the dof vector type */
	typedef Eigen::Matrix<unsigned, 1, quantity_dimension * nset_t::num_nodes> dofs_t;
	/** \brief indicates that a field view is not dirac */
	static bool const is_dirac = false;
};

/**
 * \brief Specialisation of class field_view for the isoparametric field view case
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType, class Dimension>
class field_impl<field_view<ElemType, field_option::isoparametric, Dimension> > :
	public ElemType
{
private:
	/** \brief the degree of freedom vector type */
	typedef field_traits<field_view<ElemType, field_option::isoparametric, Dimension> > traits_t;

public:
	/**
	 * \brief return underlying element
	 * \return the element of the field
	 */
	ElemType const &get_elem(void) const
	{
		// static cast
		return *this;
	}

	/** \brief return DOF vector
	 * TODO RETURN REFERENCE OF LOCAL HERE
	 */
	typename traits_t::dofs_t const &get_dofs(void) const
	{
		return ElemType::get_nodes();
	}
};


/** \brief Specialisation of field_traits for the constant field case */
template <class ElemType, class Dimension>
struct field_traits<field_view<ElemType, field_option::constant, Dimension> >
{
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief type of N-set */
	typedef constant_shape_set<typename elem_t::domain_t> nset_t;

	enum { quantity_dimension = Dimension::value };
	/** \brief the dof vector type */
	typedef Eigen::Matrix<unsigned, 1, quantity_dimension * nset_t::num_nodes> dofs_t;
	/** \brief indicates that a field view is not dirac */
	static bool const is_dirac = false;
};

/**
 * \brief Specialisation of class ::field_impl for the case of a constant field view
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType, class Dimension>
class field_impl<field_view< ElemType, field_option::constant, Dimension> >:
	public ElemType
{
private:
	typedef field_traits<field_view<ElemType, field_option::constant, Dimension> > traits_t;

public:
	/** \brief return underlying element */
	ElemType const &get_elem(void) const
	{
		return *this;
	}

	/** \brief return DOF vector */
	typename traits_t::dofs_t const &get_dofs(void) const
	{
		return ElemType::get_id();
	}
};


/**
 * \brief Field automatically generated from an element using a field generation option
 * \tparam ElemType the element type the field is associated with
 * \tparam Option field_option::constant or field_option::isoparametric
 */
template <class ElemType, class Option, class Dimension>
class field_view :
	public field_base<field_view<ElemType, Option, Dimension> >,
	public field_impl<field_view<ElemType, Option, Dimension> >
{
public:
	/** \brief the crtp base type */
	typedef field_base<field_view<ElemType, Option, Dimension> > crtp_base_t;
	/** \brief the implementation class type */
	typedef field_impl<field_view<ElemType, Option, Dimension> > impl_t;

	/** \brief self-returning metafunction */
	typedef field_view type;

	/** \brief the element type shorthand */
	typedef typename crtp_base_t::elem_t elem_t;

	enum {
		id = crtp_base_t::id
	};

	using impl_t::get_elem;
	using impl_t::get_dofs;
};


/** \brief field class that stores the dof vector and a reference to the element */
template <class ElemType, class NSet, class Dimension>
class field;


/** \brief Specialisation of field_traits for the field_extension class */
template <class ElemType, class NSet, class Dimension>
struct field_traits<field<ElemType, NSet, Dimension> >
{
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief the N-set type */
	typedef NSet nset_t;
	/** \brief the dimensionality of the physical quantity */
	enum { quantity_dimension = Dimension::value };
	/** \brief the dof vector type */
	typedef Eigen::Matrix<unsigned, 1, quantity_dimension * nset_t::num_nodes> dofs_t;
	/** \brief indicates that a field is not a dirac field */
	static bool const is_dirac = false;
};


/**
 * \brief the field class that stores the dof vector and the element by value
 * \tparam ElemType the underlying element type
 * \tparam NSet the shape function set
 */
template <class ElemType, class NSet, class Dimension>
class field_impl<field<ElemType, NSet, Dimension> >
{
public:
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief the dofs vector type */
	typedef field_traits<field<ElemType, NSet, Dimension> > traits_t;
	/** \brief the dof vector type */
	typedef typename traits_t::dofs_t dofs_t;

	/** \brief constructor from element and dof vector
	 * \param [in] elem the element reference
	 * \param [in] dofs the dof vector reference
	 */
	field_impl(element_base<elem_t> const &elem, dofs_t const &dofs) :
		m_elem(elem.derived()), m_dofs(dofs)
	{
	}

	/** \brief return underlying element */
	elem_t const &get_elem(void) const
	{
		return m_elem;
	}

	/** \return dofs */
	dofs_t const &get_dofs(void) const
	{
		return m_dofs;
	}

protected:
	/** \brief the element part by value */
	elem_t m_elem;
	/** \brief the dofs vector part by value */
	dofs_t m_dofs;
};

/**
 * \brief the field class that stores the dof vector and the element by value
 * \tparam ElemType the underlying element type
 * \tparam NSet the shape function set
 */
template <class ElemType, class NSet, class Dimension>
class field :
	public field_base<field<ElemType, NSet, Dimension> >,
	public field_impl<field<ElemType, NSet, Dimension> >
{
public:
	/** \brief self-returning metafunction */
	typedef field type;
	/** \brief the CRTP base type */
	typedef field_base<field<ElemType, NSet, Dimension> > base_t;
	/** \brief the implementation class type */
	typedef field_impl<field<ElemType, NSet, Dimension> > impl_t;

	/** \brief the element type */
	typedef typename base_t::elem_t elem_t;
	/** \brief the dofs vector type */
	typedef typename base_t::dofs_t dofs_t;

	using impl_t::get_elem;
	using impl_t::get_dofs;

	/** \brief constructor from element and dof vector
	 * \param [in] elem the element reference
	 * \param [in] dofs the degree of freedom vector reference
	 */
	field(element_base<elem_t> const &elem, dofs_t const &dofs) :
		impl_t(elem, dofs)
	{
	}
};

/** \brief field view factory */
template <class Elem, class Option, class Dimension>
field_view<Elem, Option, Dimension> const &
 create_field_view(element_base<Elem> const & e, Option, Dimension)
{
	return static_cast<field_view<Elem, Option, Dimension> const &>(e.derived());
}


/** \brief constant field view factory */
template <class Elem, class Dimension>
field_view<Elem, field_option::constant, Dimension> const &
	constant_view(element_base<Elem> const & e, Dimension)
{
	return create_field_view(e, field_option::constant(), Dimension());
}

/** \brief isoparametric field view factory */
template <class Elem, class Dimension>
field_view<Elem, field_option::isoparametric, Dimension> const &
	isoparametric_view(element_base<Elem> const & e, Dimension)
{
	return create_field_view(e, field_option::isoparametric(), Dimension());
}

/** \brief dirac field view factory */
template <class Field>
dirac_field<Field> const &
	dirac(field_base<Field> const & f)
{
	return static_cast<dirac_field<Field> const &>(
		static_cast<field_impl<Field> const &>(f.derived()));
}

#endif

