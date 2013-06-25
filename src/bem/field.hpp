/**
 * \file field.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class ::field and its specialisations
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"

struct isoparametric_field;	///< \brief tag to describe an isoparametric field
struct constant_field;		///< \brief tag to describe a constant field

template <class Derived>
struct field_traits;

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
	/** \brief the dofs vector type */
	typedef typename traits_t::dofs_t dofs_t;

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

	/**
	 * \brief return dofs vector
	 * \return the vector containing the degrees of freedoms
	 */
	dofs_t const &get_dofs(void) const
	{
		return derived().get_dofs();
	}

protected:
	elem_t const &m_elem;	///< \brief the element the field is attached to
};


/**
 * \brief A field is an element extended with a shape function set.
 * \details In our present implementation, the shape function set is not defined generally,
 * but is determined by the template parameter FieldOption: constant_field or
 * isoparametric_field. with this simplification, a field is just a proxy that
 * refers to an element.
 * \tparam ElemType the element type the field is associated with
 * \tparam FieldOption constant_field or isoparametric_field
 */
template <class ElemType, class FieldOption>
class field_view;

/** \brief Specialisation of field_traits for the isoparametric field case */
template <class ElemType>
struct field_traits<field_view<ElemType, isoparametric_field> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	typedef typename elem_t::nodes_t dofs_t;	/**< \brief the dof vector type */
};

/**
 * \brief Specialisation of class Field for the case of an isoparameteric field
 * \details On an isoparametric field the shape function set equals the geometrical L-set.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_view<ElemType, isoparametric_field>
	: public field_base<field_view<ElemType, isoparametric_field> >
{
public:
	/** \brief base's type */
	typedef field_base<field_view<ElemType, isoparametric_field> > base_t;

	/// \brief the field's elem type
	typedef typename base_t::elem_t elem_t;
	/// \brief the degree of freedom vector type
	typedef typename base_t::dofs_t dofs_t;

	typedef isoparametric_field field_option_t;	///< \brief the field generating option type

	typedef typename elem_t::lset_t nset_t;				///< \brief N-set = L-set
	static unsigned const num_dofs = nset_t::num_nodes; ///< \brief the number of dofs

	/**
	 * \brief constructor passing argument to base constructor
	 * \param [in] elem constant reference to underlying element
	 */
	field_view(elem_t const &elem) : base_t(elem)
	{
	}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return base_t::get_elem().get_nodes();
	}
};


/** \brief Specialisation of field_traits for the constant field case */
template <class ElemType>
struct field_traits<field_view<ElemType, constant_field> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	typedef typename elem_t::id_t dofs_t;	/**< \brief the dof vector type */
};

/**
 * \brief Specialisation of class Field for the case of a constant field
 * \details On a constant field the shape function set equals the constant shape
 * function set of the domain.
 * \tparam ElemType the element type the field is associated with
 */
template <class ElemType>
class field_view<ElemType, constant_field>
	: public field_base<field_view<ElemType, constant_field> >
{
public:
	/** \brief base's type */
	typedef field_base<field_view<ElemType, constant_field> > base_t;

	typedef constant_field field_option_t; ///< \brief template argument as nested type

	/// \brief the field's elem type
	typedef typename base_t::elem_t elem_t;
	/// \brief the dof vector type
	typedef typename base_t::dofs_t dofs_t;

	typedef constant_shape_set<typename elem_t::domain_t> nset_t; ///< \brief type of N-set
	static unsigned const num_dofs = nset_t::num_nodes;	///< \brief the number of dofs

	/**
	 * \brief constructor passing argument to base constructor
	 * \param [in] elem constant reference to underlying element
	 */
	field_view(elem_t const &elem) : base_t(elem)
	{
	}

	/**
	 * \brief return DOF vector
	 * \return DOF vector
	 */
	dofs_t const &get_dofs(void) const
	{
		return base_t::get_elem().get_id();
	}
};




template <class ElemType, class NSet>
class field;

/** \brief Specialisation of field_traits for the field case */
template <class ElemType, class NSet>
struct field_traits<field<ElemType, NSet> >
{
	typedef ElemType elem_t;	/**< \brief the element type */
	typedef Eigen::Matrix<unsigned, 1, NSet::num_nodes> dofs_t;	/**< \brief the dof vector type */
};


template <class ElemType, class NSet>
class field : public field_base<field<ElemType, NSet> >
{
public:
	typedef field_base<field<ElemType, NSet> > base_t;

	typedef typename base_t::elem_t elem_t;
	typedef typename base_t::dofs_t dofs_t;


	field(elem_t const &elem, dofs_t const &dofs)
		: base_t(elem), m_dofs(dofs)
	{
	}

	dofs_t const &get_dofs(void) const
	{
		return m_dofs;
	}

protected:
	dofs_t m_dofs;
};


#endif

