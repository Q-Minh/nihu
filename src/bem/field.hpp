/**
 * \file field.hpp
 * \ingroup funcspace
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class ::field and its specialisations
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

/** \brief traits class of all fields */
template <class Derived>
struct field_traits;

/** \brief metafunction assigning an id to a field
* \tparam field_t the field type
*/
template <class field_t>
struct field_id
{
	/** \brief the default field id */
	static unsigned const value =
		elem_id<typename field_traits<field_t>::elem_t>::value * 100 +
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

	/** \brief the traits class */
	typedef field_traits<Derived> traits_t;
	/** \brief the element type */
	typedef typename traits_t::elem_t elem_t;
	/** \brief the nset type */
	typedef typename traits_t::nset_t nset_t;
	/** \brief the dofs vector type */
	typedef typename traits_t::dofs_t dofs_t;
	/** \brief the number of dofs */
	static unsigned const num_dofs = nset_t::num_nodes;

	/** \brief the field identifier */
	static unsigned const id;
};

/** \brief definition of the default field identifier */
template <class Derived>
unsigned const field_base<Derived>::id = field_id<Derived>::value;


/** \brief imlementation class of a general field */
template <class Derived>
class field_impl;


template <class Field>
class dirac_field;

template <class Field>
struct field_traits<dirac_field<Field> >
{
	typedef typename field_traits<Field>::elem_t elem_t;
	typedef typename field_traits<Field>::nset_t nset_t;
	typedef typename field_traits<Field>::dofs_t dofs_t;
	static bool const is_dirac = true;
};


template <class Field>
class dirac_field :
	public field_base<dirac_field<Field> >,
	public field_impl<Field>
{
public:
	typedef dirac_field type;
	typedef typename field_impl<Field>::elem_t elem_t;
	typedef Field field_t;
};


/**
 * \brief Field automatically generated from an element using a field generation option
 * \details FieldOption can be field_option::constant or field_option::isoparametric. With this simplification,
 * a field is just a proxy that refers to an element.
 * \tparam ElemType the element type the field is associated with
 * \tparam FieldOption field_option::constant or field_option::isoparametric
 */
template <class ElemType, class FieldOption>
class field_view;

/** \brief Specialisation of field_traits for the isoparametric field view case */
template <class ElemType>
struct field_traits<field_view<ElemType, field_option::isoparametric> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	typedef typename elem_t::lset_t nset_t;	/**< \brief the dof vector type */
	typedef typename elem_t::nodes_t dofs_t;	/**< \brief the dof vector type */
	static bool const is_dirac = false;
};

/**
 * \brief Specialisation of class field_view for the isoparametric field view case
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_impl<field_view<ElemType, field_option::isoparametric> > :
	public ElemType
{
public:
	/** \brief the field's elem type */
	typedef ElemType elem_t;
	/** \brief the degree of freedom vector type */
	typedef field_traits<field_view<ElemType, field_option::isoparametric> > traits_t;
	typedef typename traits_t::dofs_t dofs_t;

	/**
	 * \brief return underlying element
	 * \return the element of the field
	 */
	elem_t const &get_elem(void) const
	{
		// static cast
		return *this;
	}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	typename traits_t::dofs_t const &get_dofs(void) const
	{
		return elem_t::get_nodes();
	}
};


/** \brief Specialisation of field_traits for the constant field case */
template <class ElemType>
struct field_traits<field_view<ElemType, field_option::constant> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	/** \brief type of N-set */
	typedef constant_shape_set<typename elem_t::domain_t> nset_t; 
	typedef typename elem_t::id_t dofs_t;	/**< \brief the dof vector type */
	static bool const is_dirac = false;
};

/**
 * \brief Specialisation of class field_impl for the case of a constant field view
 * \details On a constant field the shape function set equals the constant shape
 * function set of the domain.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_impl<field_view< ElemType, field_option::constant> >:
	public ElemType
{
public:
	/** \brief the field's elem type */
	typedef ElemType elem_t;
	/** \brief the degree of freedom vector type */
	typedef field_traits<field_view<ElemType, field_option::constant> > traits_t;
	typedef typename traits_t::dofs_t dofs_t;

	/**
	 * \brief return underlying element
	 * \return the element of the field
	 */
	elem_t const &get_elem(void) const
	{
		return *this;
	}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return elem_t::get_id();
	}
};


/**
 * \brief Specialisation of class field_view for the isoparametric field view case
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType, class Option>
class field_view :
	public field_base<field_view<ElemType, Option> >,
	public field_impl<field_view<ElemType, Option> >
{
public:
	typedef field_view type;
	typedef typename field_base<field_view<ElemType, Option> >::elem_t elem_t;
};


/** \brief field class that stores the dof vector and a reference to the element */
template <class ElemType, class NSet>
class field;


/** \brief Specialisation of field_traits for the field_extension class */
template <class ElemType, class NSet>
struct field_traits<field<ElemType, NSet> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	typedef NSet nset_t;
	typedef Eigen::Matrix<unsigned, 1, nset_t::num_nodes> dofs_t;	/**< \brief the dof vector type */
	static bool const is_dirac = false;
};


/**
* \brief the field class that stores the dof vector and the element by value
* \tparam ElemType the underlying element type
* \tparam NSet the shape function set

*/
template <class ElemType, class NSet>
class field_impl<field<ElemType, NSet> > 
{
public:
	/** \brief the element type */
	typedef ElemType elem_t;
	/** \brief the dofs vector type */
	typedef field_traits<field<ElemType, NSet> > traits_t;

	typedef typename traits_t::dofs_t dofs_t;

	field_impl(element_base<elem_t> const &elem, dofs_t const &dofs) :
		m_elem(elem.derived()), m_dofs(dofs)
	{
	}

	/**
	 * \brief return underlying element
	 * \return the element of the field
	 */
	elem_t const &get_elem(void) const
	{
		return m_elem;
	}

	/**
  	 * \return dofs
 	 */
	dofs_t const &get_dofs(void) const
	{
		return m_dofs;
	}

protected:
	elem_t m_elem;
	dofs_t m_dofs;
};

/**
* \brief the field class that stores the dof vector and the element by value
* \tparam ElemType the underlying element type
* \tparam NSet the shape function set
*/
template <class ElemType, class NSet>
class field : 
	public field_base<field<ElemType, NSet> >,
	public field_impl<field<ElemType, NSet> >
{
public:
	/** \brief the CRTP base type */
	typedef field_base<field<ElemType, NSet> > base_t;
	typedef field_impl<field<ElemType, NSet> > impl_t;

	/** \brief the element type */
	typedef typename base_t::elem_t elem_t;
	/** \brief the dofs vector type */
	typedef typename base_t::dofs_t dofs_t;


	field(element_base<elem_t> const &elem, dofs_t const &dofs) :
		impl_t(elem, dofs)
	{
	}
};

/**
 * \brief field view factory
 */
template <class Elem, class Option> 
field_view<Elem, Option> const &
 create_field_view(element_base<Elem> const & e, Option)
{
	return static_cast<field_view<Elem, Option> const &>(e.derived());
}


/**
 * \brief constant field view factory
 */
template <class Elem>
field_view<Elem, field_option::constant> const &
	constant_view(element_base<Elem> const & e)
{
	return create_field_view(e, field_option::constant());
}

/**
 * \brief isoparametric field view factory
 */
template <class Elem>
field_view<Elem, field_option::isoparametric> const &
	isoparametric_view(element_base<Elem> const & e)
{
	return create_field_view(e, field_option::isoparametric());
}

/**
 * \brief dirac field view factory 
 */
template <class Field>
dirac_field<Field> const & 
	dirac(field_base<Field> const & f)
{
	return static_cast<dirac_field<Field> const &>(
		static_cast<field_impl<Field> const &>(f.derived()));
}

#endif

