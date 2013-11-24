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
 * \file field_type_accelerator.hpp
 * \ingroup quadrature
 * \brief Acceleration of field types by storing shape functions for each quadrature point
 */

#ifndef FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

#include "global_definitions.hpp"
#include "field_type_acceleration_option.hpp"
#include "field.hpp"
#include "quadrature.hpp"
#include "../util/casted_iterator.hpp"
#include "../util/dual_range.hpp"
#include "../util/store_pattern.hpp"
#include "../util/pool_pattern.hpp"


/** \brief Stores a quadrature point and a shape function vector
 * \tparam Field the field type to accelerate
 * \tparam Family the quadrature family
 * \tparam Acceleration hard or soft
 */
template<class Field, class Family, class Acceleration>
class field_type_accelerator_elem;

/** \brief specialisation of ::field_type_accelerator_elem for the soft case */
template <class Field, class Family>
class field_type_accelerator_elem<Field, Family, acceleration::soft> :
	public quadrature_elem<
		typename Field::elem_t::domain_t::xi_t,
		typename Field::elem_t::domain_t::scalar_t
	>
{
public:
	/** \brief the base type */
	typedef quadrature_elem<
		typename Field::elem_t::domain_t::xi_t, typename Field::domain_t::scalar_t
	> base_t;

	/** \brief return shape function
	 * \return the shape function vector
	 */
	typename Field::nset_t::shape_t get_N(void) const
	{
		return Field::nset_t::eval_shape(base_t::get_xi());
	}
};


/** \brief specialisation of ::field_type_accelerator_elem for the hard acceleration case */
template <class Field, class Family>
class field_type_accelerator_elem<Field, Family, acceleration::hard> :
	public quadrature_elem<
		typename Field::elem_t::domain_t::xi_t,
		typename Field::elem_t::domain_t::scalar_t
	>
{
public:
	/** \brief the base type */
	typedef quadrature_elem<
		typename Field::elem_t::domain_t::xi_t,
		typename Field::elem_t::domain_t::scalar_t
	> base_t;

	/** \brief constructor computing the shape function
	 * \param [in] base the quadrature element reference
	 */
	field_type_accelerator_elem(base_t const &base) :
		base_t(base),
		m_nset(Field::nset_t::eval_shape(base_t::get_xi()))
	{
	}

	/** \brief return shape function
	 * \return the shape function vector
	 */
	typename Field::nset_t::shape_t const &get_N(void) const
	{
		return m_nset;
	}

protected:
	/** \brief the accelerated n-set */
	typename Field::nset_t::shape_t m_nset;
};



// forward declaration
template <class Field, class Family, class Acceleration, class Enable = void>
class field_type_accelerator;

/** \brief field type accelerator for the hard non-dirac case
 * \tparam Field the field type to accelerate
 * \tparam the quadrature family type
 */
template <class Field, class Family>
class field_type_accelerator<
	Field, Family, acceleration::hard,
	typename std::enable_if<!field_traits<Field>::is_dirac>::type
	> :
	public EigenStdVector<
		field_type_accelerator_elem<Field, Family, acceleration::hard>
	>::type
{
public:
	/** \brief the accelerator element type */
	typedef field_type_accelerator_elem<
		Field, Family, acceleration::hard
	> accelerator_elem_t;
	/** \brief the base type */
	typedef typename EigenStdVector<accelerator_elem_t>::type base_t;

	/** \brief the quadrature type */
	typedef typename quadrature_type<
		Family, typename Field::elem_t::domain_t
	>::type quadrature_t;

	/** \brief contructor from a quadrature
	 * \param [in] quadrature the quadrature instance
	 */
	field_type_accelerator(quadrature_t const &quadrature)
	{
		this->reserve(quadrature.size());
		for (unsigned i = 0; i < quadrature.size(); ++i)
			this->push_back(accelerator_elem_t(quadrature[i]));
	}

	/** \brief constructor from quadrature order
	 * \param [in] order the quadrature order
	 */
	field_type_accelerator(unsigned order) :
		field_type_accelerator(quadrature_t(order))
	{
	}
};


/** \brief field type accelerator for the soft non-dirac case
 * \tparam Field the field type to accelerate
 * \tparam the quadrature family type
 */
template <class Field, class Family>
class field_type_accelerator<Field, Family, acceleration::soft,
	typename std::enable_if<!field_traits<Field>::is_dirac>::type
	> :
	public quadrature_type<Family, typename Field::elem_t::domain_t>::type
{
public:
	/** \brief the base type */
	typedef typename quadrature_type<Family, typename Field::elem_t::domain_t>::type base_t;
	/** \brief the accelerator element type */
	typedef field_type_accelerator_elem<Field, Family, acceleration::soft> accelerator_elem_t;
	/** \brief the iterator type */
	typedef casted_iterator<
		typename base_t::const_iterator,
		accelerator_elem_t
	> const_iterator;

	/** \brief return begin iterator
	 * \return begin iterator
	 */
	const_iterator begin(void) const
	{
		return base_t::cbegin();
	}

	/** \brief return end iterator
	 * \return end iterator
	 */
	const_iterator end(void) const
	{
		return base_t::cend();
	}
};


/** \brief an index class defined to use as a base class */
struct index_t
{
	/** \brief constructor */
	index_t(unsigned idx) :	m_idx(idx) {}
	/** \brief the underlying index value */
	unsigned m_idx;
};


/** \brief accelerator elem for a Dirac field
 * \tparam NSet the shape functin set type
 */
template <class NSet>
class dirac_field_type_accelerator_elem :
	public index_t
{
public:
	/** \brief return quadrature weight
	 * \return quadrature weight (constant 1.0)
	 */
	constexpr typename NSet::scalar_t get_w(void) const
	{
		return 1.0;
	}

	/** \brief return shape set vector
	 * \return constant shape set vector
	 */
	constexpr typename NSet::shape_t get_N(void) const
	{
		return NSet::shape_t::Unit(index_t::m_idx);
	}

	/** \brief return quadrature point location
	 * \return constant quadrature point location
	 */
	constexpr typename NSet::xi_t get_xi(void) const
	{
		return NSet::corner_at(index_t::m_idx);
	}
};

/** \brief an iterator class used for the dirac accelerators */
class dirac_field_type_accelerator_iterator
{
public:
	/** \brief constructor
	 * \param [in] idx the index
	 */
	dirac_field_type_accelerator_iterator(index_t const &idx) :
		m_idx(idx)
	{
	}

	/** \brief increment operator
	 * \return the incremented iterator reference
	 */
	dirac_field_type_accelerator_iterator &operator++(void)
	{
		++m_idx.m_idx;
		return *this;
	}

	/** \brief check inequality
	 * \param [in] other the other iterator
	 * \return true if the two iterators are different
	 */
	bool operator !=(dirac_field_type_accelerator_iterator const &other) const
	{
		return m_idx.m_idx != other.m_idx.m_idx;
	}

	/** \brief check equality
	 * \param [in] other the other iterator
	 * \return true if the two are equal
	 */
	bool operator ==(dirac_field_type_accelerator_iterator const &other) const
	{
		return !(*this != other);
	}

	/** \brief dereference operator
	 * \return the dereferred index
	 */
	index_t const &operator*(void) const
	{
		return m_idx;
	}

	/** \brief pointer dereference operator
	 * \return the address of the stored index
	 */
	index_t const *operator->(void) const
	{
		return &m_idx;
	}

private:
	index_t m_idx;
};


/** \brief specialisation of ::field_type_accelerator for the Dirac field type
 * \tparam Field the field type
 * \tparam Family the quadrature family
 * \tparam Acceleration the acceleration option
 */
template <class Field, class Family, class Acceleration>
class field_type_accelerator<Field, Family, Acceleration,
	typename std::enable_if<field_traits<Field>::is_dirac>::type
	>
{
public:
	typedef dirac_field_type_accelerator_elem<
		typename Field::nset_t
	> accelerator_elem_t;

	/** \brief default constructor */
	field_type_accelerator(unsigned)
	{
	}

	typedef casted_iterator<
		dirac_field_type_accelerator_iterator,
		accelerator_elem_t
	> const_iterator;

	/** \brief return begin iterator
	 * \return begin iterator
	 */
	constexpr static const_iterator begin(void)
	{
		return dirac_field_type_accelerator_iterator(index_t(0));
	}

	/** \brief return end iterator
	 * \return end iterator
	 */
	constexpr static const_iterator end(void)
	{
		return dirac_field_type_accelerator_iterator(index_t(Field::num_dofs));
	}
};


template <class Field, class Family, class Acceleration, unsigned MaxOrder>
class field_type_accelerator_pool :
	public pool<field_type_accelerator<Field, Family, acceleration::hard>, MaxOrder>
{
public:
	typedef pool<field_type_accelerator<Field, Family, acceleration::hard>, MaxOrder> base_t;

	/** \brief this user-specified default constructor is needed by clang */
	field_type_accelerator_pool() :
		base_t()
	{
	}
};


template <class Field, class Family, unsigned MaxOrder>
class field_type_accelerator_pool<Field, Family, acceleration::soft, MaxOrder> :
	public pool<typename quadrature_type<Family, typename Field::elem_t::domain_t>::type, MaxOrder>
{
public:
	typedef pool<typename quadrature_type<Family, typename Field::elem_t::domain_t>::type, MaxOrder> base_t;
	typedef field_type_accelerator<Field, Family, acceleration::soft> accelerator_t;

	/** \brief this user-specified default constructor is needed by clang */
	field_type_accelerator_pool() :
		base_t()
	{
	}

	accelerator_t const &operator[](unsigned idx) const
	{
		return static_cast<accelerator_t const &>(base_t::operator[](idx));
	}
};


template <class TestAccelerator, class TrialAccelerator, class IterationMode>
class dual_field_type_accelerator :
	public dual_range<
		IterationMode,
		typename TestAccelerator::const_iterator,
		typename TrialAccelerator::const_iterator
	>
{
public:
	typedef dual_range<
		IterationMode,
		typename TestAccelerator::const_iterator,
		typename TrialAccelerator::const_iterator
	> base_t;

	dual_field_type_accelerator(
		TestAccelerator const &test_accelerator,
		TrialAccelerator const &trial_accelerator
	) : base_t(
		test_accelerator.begin(), test_accelerator.end(),
		trial_accelerator.begin(), trial_accelerator.end())
	{
	}
};

template <class TestAccelerator, class TrialAccelerator, class IterationMode>
dual_field_type_accelerator<TestAccelerator, TrialAccelerator, IterationMode>
	create_dual_field_type_accelerator(
		TestAccelerator const &test_acc,
		TrialAccelerator const &trial_acc,
		IterationMode
	)
{
	return dual_field_type_accelerator<TestAccelerator, TrialAccelerator, IterationMode>(
		test_acc, trial_acc);
}


#endif // FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

