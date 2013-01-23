/**
 * \file field.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class Field and its specialisation
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"
#include "../tmp/general.hpp"

/** \brief tag class used to describe an isoparametric field */
struct isoparametric_field;
/** \brief tag class used to describe a constant field */
struct constant_field;

/**
 * \brief A field is an element extended with a shape function set.
 * \details In our present implementation, the shape function set is not defined generally,
 * but is determined by the template parameter FieldOption: constant_field or isoparametric_field
 * with this simplification, a field is just a proxy that refers to an element 
 * \tparam ElemType the element type the field is associated with
 * \tparam FieldOption constant_field or isoparametric_field
 */
template <class ElemType, class FieldOption>
class field;


/**
 * \brief Specialisation of class Field for the case of an isoparameteric field
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 * \tparam FieldOption constant_field or isoparametric_field
 */
template <class ElemType>
class field<ElemType, isoparametric_field>
{
public:
	/** \brief template parameter as nested type */
	typedef ElemType elem_t;
	/** \brief template parameter as nested type */
	typedef isoparametric_field field_option_t;

	/** \brief the L-set's type inherited from the element template */
	typedef typename elem_t::lset_t nset_t;
	/** \brief the number of dofs */
	static unsigned const num_dofs = nset_t::num_nodes;
	/** \brief type of the field's DOF vector */
	typedef typename elem_t::nodes_t dofs_t;

	/**
	 * \brief constructor initialising the reference member
	 * \param elem The element the field is attached to
	 */
	field(elem_t const &elem) : elem(elem)
	{
	}
	
	elem_t const &get_elem(void) const
	{
		return elem;
	}

	/**
	 * \brief return DOF vector
	 * \returns DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return elem.get_nodes();
	}

protected:
	/** \brief the element the field is attached to */
	elem_t const &elem;
};


template <class ElemType>
class field<ElemType, constant_field>
{
public:
	/** \brief template parameter as nested type */
	typedef ElemType elem_t;
	/** \brief template parameter as nested type */
	typedef constant_field field_option_t;

	/** \brief the element domain's constant shape set type is the N-set */
	typedef constant_shape_set<typename elem_t::domain_t> nset_t;
	/** \brief the number of dofs */
	static unsigned const num_dofs = nset_t::num_nodes;
	/** \brief type of the field's DOF vector */
	typedef typename elem_t::id_t dofs_t;

	/**
	 * \brief constructor initialising the reference member
	 * \param elem The element the field is attached to
	 */
	field(elem_t const &elem) : elem(elem)
	{
	}

	elem_t const &get_elem(void) const
	{
		return elem;
	}

	/**
	 * \brief return DOF vector
	 * \returns DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return elem.get_id();
	}

protected:
	/** \brief the element the field is attached to */
	elem_t const &elem;
};

#endif

