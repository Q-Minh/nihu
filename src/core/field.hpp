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

namespace NiHu
{

/** \brief collect options used to convert an element into a field view */
namespace field_option
{
	/** \brief tag to describe an isoparametric field */
	struct isoparametric {};
	/** \brief tag to describe a constant field */
	struct constant {};
}

/** \brief type indicating a 1D valued function space */
typedef std::integral_constant<unsigned, 1> _1d;
/** \brief type indicating a 2D valued function space */
typedef std::integral_constant<unsigned, 2> _2d;
/** \brief type indicating a 3D valued function space */
typedef std::integral_constant<unsigned, 3> _3d;

/** \brief assign a field to a tag
 * \tparam field_tag the tag
 */
template <class field_tag>
struct tag2field;


namespace field_traits
{
    /** \brief assigns the element type to the field */
	template <class Derived>
	struct elem_type;

    /** \brief assigns the N-set type to the field */
	template <class Derived>
	struct nset_type;

    /** \brief assigns the dimensionality of the interpolated physical quantity */
	template <class Derived>
	struct quantity_dimension;

    /** \brief indicate if the field stores its DOF vector or computes it on the fly */
	template <class Derived>
	struct is_dof_vector_stored : std::false_type {};

    /** \brief indicate if the field is a Dirac field or not */
	template <class Derived>
	struct is_dirac : std::false_type {};

	/** \brief assign a numeric ID to the field */
	template <class Derived>
	struct id
	{
		/** \brief the default field id */
		enum { value =
			elem_type<Derived>::type::id * 100 + nset_type<Derived>::type::num_nodes
		};
	};

    /** \brief assign the DOF vector value type to the field type */
	template <class Derived>
	struct dof_vector_type
	{
		typedef Eigen::Matrix<
			unsigned,
			quantity_dimension<Derived>::value * nset_type<Derived>::type::num_nodes,
			1
		> type;
	};

    /** \brief assign the DOF vector return type to the field type */
	template <class Derived>
	struct dof_vector_return_type : std::conditional<
		is_dof_vector_stored<Derived>::value,
		typename dof_vector_type<Derived>::type const &,
		typename dof_vector_type<Derived>::type
	> {};
}


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

	/** \brief the element type */
	typedef typename field_traits::elem_type<Derived>::type elem_t;
	/** \brief the nset type */
	typedef typename field_traits::nset_type<Derived>::type nset_t;
	/** \brief the dofs vector type */
	typedef typename field_traits::dof_vector_type<Derived>::type dofs_t;
	/** \brief the dofs vector return type */
	typedef typename field_traits::dof_vector_return_type<Derived>::type dofs_return_t;

	enum {
		/** \brief the number of dofs */
		num_dofs = nset_t::num_nodes,
		/** \brief the field id */
		id = field_traits::id<Derived>::value,
		/** \brief the quantity's dimensionality */
		quantity_dimension = field_traits::quantity_dimension<Derived>::value
	};

	/** \brief return underlying element */
	elem_t const &get_elem(void) const
	{
		return derived().get_elem();
	}

	/** \brief return DOF vector */
	dofs_return_t get_dofs(void) const
	{
		return derived().get_dofs();
	}
};



/** \brief implementation class of a general field */
template <class Derived>
class field_impl;


// forward declaration
template <class Field>
class dirac_field;

namespace field_traits
{
	template <class Derived>
	struct elem_type<dirac_field<Derived> > : elem_type<Derived> {};

	template <class Derived>
	struct nset_type<dirac_field<Derived> > : nset_type<Derived> {};

	template <class Derived>
	struct quantity_dimension<dirac_field<Derived> > : quantity_dimension<Derived> {};

	template <class Derived>
	struct is_dirac<dirac_field<Derived> > : std::true_type {};
}


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
	typedef typename field_traits::elem_type<dirac_field<Field> >::type elem_t;
	/** \brief store the original field type for the iterator */
	typedef Field field_t;

	using impl_t::get_elem;
	using impl_t::get_dofs;
};


// forward declaration
template <class ElemType, class FieldOption, class Dimension = _1d>
class field_view;

namespace field_traits
{
    /** \brief assign an element type to a field view */
	template <class ElemType, class FieldOption, class Dimension>
	struct elem_type<field_view<ElemType, FieldOption, Dimension> > : ElemType {};

    /** \brief assign an N-set type to a field view */
	template <class ElemType, class Dimension>
	struct nset_type<field_view<ElemType, field_option::isoparametric, Dimension> >
	{
		typedef  typename ElemType::lset_t type;
	};

    /** \brief assign an N-set type to a constant field view */
	template <class ElemType, class Dimension>
	struct nset_type<field_view<ElemType, field_option::constant, Dimension> >
	{
		typedef constant_shape_set<typename ElemType::domain_t> type;
	};

    /** \brief assign the dimension of the interpolated quantity to a field view */
	template <class ElemType, class FieldOption, class Dimension>
	struct quantity_dimension<field_view<ElemType, FieldOption, Dimension> > : Dimension {};
}


/**
 * \brief Specialisation of class field_impl for the isoparametric field view case
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType, class Dimension>
class field_impl<field_view<ElemType, field_option::isoparametric, Dimension> > :
	public ElemType
{
	typedef field_view<ElemType, field_option::isoparametric, Dimension> Derived;
	typedef typename field_traits::dof_vector_type<Derived>::type dofs_t;
	typedef typename field_traits::dof_vector_return_type<Derived>::type dofs_return_t;
public:
	/** \brief return underlying element */
	ElemType const &get_elem(void) const
	{
		return *this;	// static cast
	}

	/** \brief return DOF vector */
	dofs_return_t get_dofs(void) const
	{
	enum { D = Dimension::value };
		dofs_t dofs;
		for (unsigned n = 0; n < ElemType::num_nodes; ++n)
			for (unsigned d = 0; d < D; ++d)
				dofs(n*D+d) = ElemType::get_nodes()(n)*D + d;
		return dofs;
	}
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
	typedef field_view<ElemType, field_option::constant, Dimension> Derived;
	typedef typename field_traits::dof_vector_type<Derived>::type dofs_t;
	typedef typename field_traits::dof_vector_return_type<Derived>::type dofs_return_t;

public:
	/** \brief return underlying element */
	ElemType const &get_elem(void) const
	{
		return *this; // static cast
	}

	/** \brief return DOF vector */
	dofs_return_t get_dofs(void) const
	{
		enum { D = Dimension::value };
		return dofs_t::LinSpaced(D, ElemType::get_id()(0)*D, ElemType::get_id()(0)*D+D-1);
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
	/** \brief the crtp base type */
	typedef field_base<field_view<ElemType, Option, Dimension> > crtp_base_t;
	/** \brief the implementation class type */
	typedef field_impl<field_view<ElemType, Option, Dimension> > impl_t;

public:
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


namespace field_traits
{
    /** \brief assign an element type to a field */
	template <class ElemType, class NSet, class Dimension>
	struct elem_type<field<ElemType, NSet, Dimension> > : ElemType {};

    /** \brief assign an N-set type to a field */
	template <class ElemType, class NSet, class Dimension>
	struct nset_type<field<ElemType, NSet, Dimension> >
	{
		typedef NSet type;
	};

    /** \brief assign the dimension of the interpolated quantity to a field */
	template <class ElemType, class NSet, class Dimension>
	struct quantity_dimension<field<ElemType, NSet, Dimension> > : Dimension {};
}


/**
 * \brief the field class that stores the dof vector and the element by value
 * \tparam ElemType the underlying element type
 * \tparam NSet the shape function set
 */
template <class ElemType, class NSet, class Dimension>
class field_impl<field<ElemType, NSet, Dimension> >
{
public:
    /** \brief the derived field type */
	typedef field<ElemType, NSet, Dimension> Derived;
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief the dof vector type */
	typedef typename field_traits::dof_vector_type<Derived>::type dofs_t;
	/** \brief the dof vector type */
	typedef typename field_traits::dof_vector_return_type<Derived>::type dofs_return_t;

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
	dofs_return_t get_dofs(void) const
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
template <class ElemType, class NSet, class Dimension = _1d>
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
template <class Elem, class Option, class Dimension = _1d>
field_view<Elem, Option, Dimension> const &
	create_field_view(element_base<Elem> const & e, Option, Dimension dim = Dimension())
{
	return static_cast<field_view<Elem, Option, Dimension> const &>(e.derived());
}


/** \brief constant field view factory */
template <class Elem, class Dimension = _1d>
field_view<Elem, field_option::constant, Dimension> const &
	constant_view(element_base<Elem> const & e, Dimension dim = Dimension())
{
	return create_field_view(e, field_option::constant(), dim);
}

/** \brief isoparametric field view factory */
template <class Elem, class Dimension = _1d>
field_view<Elem, field_option::isoparametric, Dimension> const &
	isoparametric_view(element_base<Elem> const & e, Dimension dim = Dimension())
{
	return create_field_view(e, field_option::isoparametric(), dim);
}

/** \brief dirac field view factory */
template <class Field>
dirac_field<Field> const &
	dirac(field_base<Field> const & f)
{
	return static_cast<dirac_field<Field> const &>(
		static_cast<field_impl<Field> const &>(f.derived()));
}

}

#endif

