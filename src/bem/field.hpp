/**
 * \file field.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief declaration of class Field
 */
#ifndef FIELD_HPP_INCLUDED
#define FIELD_HPP_INCLUDED

#include "element.hpp"
#include "../tmp/general.hpp"

struct isoparametric_field;
struct constant_field;

template <class ElemType, class FieldOption>
class Field;

template <class ElemType>
class Field<ElemType, isoparametric_field>
{
public:
	typedef ElemType elem_t;
	typedef isoparametric_field field_option_t;

	typedef typename elem_t::lset_t nset_t;
	typedef typename elem_t::nodes_t dofs_t;

	Field(elem_t const &elem) : elem(elem)
	{
	}
	
	dofs_t const &get_dofs(void) const
	{
		return elem.get_nodes();
	}

protected:
	elem_t const &elem;
};


template <class ElemType>
class Field<ElemType, constant_field>
{
public:
	typedef ElemType elem_t;
	typedef constant_field field_option_t;

	typedef constant_shape_set<typename elem_t::domain_t> nset_t;
	typedef typename elem_t::id_t dofs_t;

	Field(elem_t const &elem) : elem(elem)
	{
	}

	dofs_t const &get_dofs(void) const
	{
		return elem.get_id();
	}

protected:
	elem_t const &elem;
};

#endif

