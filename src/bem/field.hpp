/**
 * \file field.hpp
 * \ingroup funcspace
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class ::field and its specialisations
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"

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
private:
	/** \brief helper function to return reference to the derived class */
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	/** \brief helper function to return reference to the derived class */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
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
	static unsigned const id; // = field_id<Derived>::value;
	

	/**
	 * \brief return underlying element
	 * \return the element of the field
	 */
	elem_t const &get_elem(void) const
	{
		return derived().get_elem();
	}

	/**
	 * \brief return dofs vector
	 * \return the vector containing the degrees of freedoms
	 */
	dofs_t const &get_dofs(void) const
	{
		return derived().get_dofs();
	}
};

template <class Derived>
unsigned const field_base<Derived>::id = field_id<Derived>::value;


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
};

/**
 * \brief Specialisation of class field_view for the isoparametric field view case
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_view<ElemType, field_option::isoparametric>
	: public field_base<field_view<ElemType, field_option::isoparametric> >
{
public:
	/** \brief base's type */
	typedef field_base<field_view<ElemType, field_option::isoparametric> > base_t;

	typedef field_view type;

	/** \brief the field's elem type */
	typedef typename base_t::elem_t elem_t;
	/** \brief the degree of freedom vector type */
	typedef typename base_t::dofs_t dofs_t;

	/**
	 * \brief constructor simply passing argument to base constructor
	 * \param [in] elem constant reference to underlying element
	 */
	field_view(elem_t const &elem) : m_elem(elem)
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
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return m_elem.get_nodes();
	}
	
private:
	elem_t const &m_elem;
};


/** \brief Specialisation of field_traits for the constant field case */
template <class ElemType>
struct field_traits<field_view<ElemType, field_option::constant> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	/** \brief type of N-set */
	typedef constant_shape_set<typename elem_t::domain_t> nset_t; 
	typedef typename elem_t::id_t dofs_t;	/**< \brief the dof vector type */
};

/**
 * \brief Specialisation of class Field for the case of a constant field
 * \details On a constant field the shape function set equals the constant shape
 * function set of the domain.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_view<ElemType, field_option::constant>
	: public field_base<field_view<ElemType, field_option::constant> >
{
public:
	/** \brief base's type */
	typedef field_base<field_view<ElemType, field_option::constant> > base_t;

	typedef field_view type;

	/** \brief the element type */
	typedef typename base_t::elem_t elem_t;
	/** \brief the dof vector type */
	typedef typename base_t::dofs_t dofs_t;

	/**
	 * \brief constructor passing argument to base constructor
	 * \param [in] elem constant reference to underlying element
	 */
	field_view(elem_t const &elem) : m_elem(elem)
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
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return m_elem.get_id();
	}
	
private:
	elem_t const &m_elem;
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
};

/**
* \brief the field class that stores the dof vector and the element by value
* \tparam ElemType the underlying element type
* \tparam NSet the shape function set
*/
template <class ElemType, class NSet>
class field : public field_base<field<ElemType, NSet> >
{
public:
	/** \brief the CRTP base type */
	typedef field_base<field<ElemType, NSet> > base_t;

	/** \brief the element type */
	typedef typename base_t::elem_t elem_t;
	/** \brief the dofs vector type */
	typedef typename base_t::dofs_t dofs_t;


	field(elem_t const &elem, dofs_t const &dofs) : m_elem(elem), m_dofs(dofs)
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

	dofs_t const &get_dofs(void) const
	{
		return m_dofs;
	}

protected:
	elem_t m_elem;
	dofs_t m_dofs;
};


#endif

