/**
 * \file field.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class Field and its specialisation
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"

/** \brief tag class used to describe an isoparametric field */
struct isoparametric_field;
/** \brief tag class used to describe a constant field */
struct constant_field;

/**
 * \brief a comon base type for both fields
 * \tparam ElemType the element type of the field
 */
template <class ElemType>
class field_base
{
public:
	typedef ElemType elem_t; /**< \brief template parameter as nested type */

	/**
	 * \brief constructor initialising the reference member
	 * \param elem The element the field is attached to
	 */
	field_base(elem_t const &elem) : m_elem(elem)
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

protected:
	/** \brief the element the field is attached to */
	elem_t const &m_elem;
};


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
 */
template <class ElemType>
class field<ElemType, isoparametric_field> : public field_base<ElemType>
{
public:
	typedef field_base<ElemType> base;	/**< \brief base's type */

	typedef typename base::elem_t elem_t;	/**< \brief the field's elem type */
	typedef isoparametric_field field_option_t;	/**< \brief the field generating option type */

	typedef typename elem_t::lset_t nset_t;				/**< \brief N-set = L-set */
	static unsigned const num_dofs = nset_t::num_nodes; /**< \brief the number of dofs */
	typedef typename elem_t::nodes_t dofs_t;			/**< \brief type of DOF vector */

	/**
	 * \brief constructor passing argument to base constructor
	 * \param elem constant reference to underlying element
	 */
	field(elem_t const &elem) : field_base<elem_t>(elem) {}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return base::m_elem.get_nodes();
	}
};


/**
 * \brief Specialisation of class Field for the case of a constant field
 * \details On a constant field the shape function set equals the constant shape function set of the domain.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field<ElemType, constant_field> : public field_base<ElemType>
{
public:
	typedef field_base<ElemType> base;	/**< \brief base's type */

	typedef typename base::elem_t elem_t;	/**< \brief the field's elem type */
	typedef constant_field field_option_t;		/**< \brief the field generating option */

	typedef constant_shape_set<typename elem_t::domain_t> nset_t; /**< \brief type of N-set */
	static unsigned const num_dofs = nset_t::num_nodes;	/**< \brief the number of dofs */
	typedef typename elem_t::id_t dofs_t;				/**< \brief type of the field's DOF vector */

	/**
	 * \brief constructor passing argument to base constructor
	 * \param elem constant reference to underlying element
	 */
	field(elem_t const &elem) : field_base<elem_t>(elem)
	{
	}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return base::m_elem.get_id();
	}
};

#endif

