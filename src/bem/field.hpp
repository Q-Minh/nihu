/**
 * \file field.hpp
 * \brief declaration of class Field
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"
#include "../tmp/general.hpp"

template <class ElemType, class NSet>
class Field
{
public:
	typedef ElemType elem_t;
	typedef NSet nset_t;

	static const unsigned num_dofs = nset_t::num_nodes;
	typedef Eigen::Matrix<unsigned, 1, num_dofs> dofs_t;

	static_assert(is_same<typename elem_t::domain_t, typename nset_t::domain_t>::value,
		"Element and NSet domains must match");

	Field() : elem(NULL)
	{
	}

	Field(elem_t const *elem, dofs_t const &dofs) : elem(elem), dofs(dofs)
	{
	}

protected:
	elem_t const *elem; // pointer needed in order to have a default constructor (std::vector)
	dofs_t dofs;
};


template <class ElemType>
class IsoParametricField
{
public:
	typedef ElemType elem_t;
	typedef typename elem_t::lset_t nset_t;

	static const unsigned num_dofs = nset_t::num_nodes;
	typedef Eigen::Matrix<unsigned, 1, num_dofs> dofs_t;

	IsoParametricField() : elem(NULL)
	{
	}

	IsoParametricField(elem_t const *elem) : elem(elem)
	{
	}

	dofs_t const &get_dofs(void) const
	{
		return elem->get_nodes();
	}

protected:
	elem_t const *elem; // pointer needed in order to have a default constructor (std::vector)
};


template <class ElemType>
class ConstantField
{
public:
	typedef ElemType elem_t;
	typedef constant_shape_set<typename elem_t::domain_t> nset_t;

	static const unsigned num_dofs = nset_t::num_nodes; // 1
	typedef unsigned dofs_t;

	ConstantField() : elem(NULL)
	{
	}

	ConstantField(elem_t const *elem) : elem(elem)
	{
	}

	dofs_t const &get_dofs(void) const
	{
		return elem->get_id();
	}

protected:
	elem_t const *elem; // pointer needed in order to have a default constructor (std::vector)
};


#endif
