/**
 * \file field_type_accelerator.hpp
 * \ingroup quadrature
 * \brief Acceleration of field types by storing shape functions for each quadrature point
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

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
		typename Field::domain_t::xi_t,
		typename Field::domain_t::scalar_t
	>
{
public:
	typedef quadrature_elem<
		typename Field::domain_t::xi_t, typename Field::domain_t::scalar_t
	> base_t;

	/** \brief return shape function
	 * \return the shape function vector
	 */
	typename Field::nset_t::shape_t get_N(void) const
	{
		return Field::nset_t::eval_shape(base_t::get_xi());
	}
};


/** \brief specialisation of ::field_type_accelerator_elem for the hard case */
template <class Field, class Family>
class field_type_accelerator_elem<Field, Family, acceleration::hard> :
	public quadrature_elem<
		typename Field::domain_t::xi_t,
		typename Field::domain_t::scalar_t
	>
{
public:
	typedef quadrature_elem<
		typename Field::domain_t::xi_t,
		typename Field::domain_t::scalar_t
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
	typename Field::nset_t::shape_t m_nset;
};



template <class Field, class Family, class Acceleration, class Enable = void>
class field_type_accelerator;

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
	typedef field_type_accelerator_elem<
		Field, Family, acceleration::hard
	> accelerator_elem_t;
	typedef typename EigenStdVector<accelerator_elem_t>::type base_t;

	typedef typename quadrature_type<
		Family, typename Field::domain_t
	>::type quadrature_t;

	field_type_accelerator(quadrature_t const &quadrature)
	{
		this->reserve(quadrature.size());
		for (unsigned i = 0; i < quadrature.size(); ++i)
			this->push_back(accelerator_elem_t(quadrature[i]));
	}

	field_type_accelerator(unsigned order) :
		field_type_accelerator(quadrature_t(order))
	{
	}
};


template <class Field, class Family>
class field_type_accelerator<Field, Family, acceleration::soft,
	typename std::enable_if<!field_traits<Field>::is_dirac>::type
	> :
	public quadrature_type<Family, typename Field::domain_t>::type
{
public:
	typedef typename quadrature_type<Family, typename Field::domain_t>::type base_t;

	typedef field_type_accelerator_elem<Field, Family, acceleration::soft> accelerator_elem_t;

	typedef casted_iterator<
		typename base_t::const_iterator,
		accelerator_elem_t
	> const_iterator;

	const_iterator begin(void) const
	{
		return base_t::cbegin();
	}

	const_iterator end(void) const
	{
		return base_t::cend();
	}
};


struct index_t
{
	index_t(unsigned idx) :	m_idx(idx) {}
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

class dirac_field_type_accelerator_iterator
{
public:
	dirac_field_type_accelerator_iterator(index_t const &idx) :
		m_idx(idx)
	{
	}

	dirac_field_type_accelerator_iterator &operator++(void)
	{
		++m_idx.m_idx;
		return *this;
	}

	bool operator !=(dirac_field_type_accelerator_iterator const &other) const
	{
		return m_idx.m_idx != other.m_idx.m_idx;
	}

	bool operator ==(dirac_field_type_accelerator_iterator const &other) const
	{
		return !(*this != other);
	}

	index_t const &operator*(void) const
	{
		return m_idx;
	}

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

	field_type_accelerator(unsigned)
	{
	}

	typedef casted_iterator<
		dirac_field_type_accelerator_iterator,
		accelerator_elem_t
	> const_iterator;

	constexpr static const_iterator begin(void)
	{
		return dirac_field_type_accelerator_iterator(index_t(0));
	}

	constexpr static const_iterator end(void)
	{
		return dirac_field_type_accelerator_iterator(index_t(Field::num_dofs));
	}
};


template <class Field, class Family, class Acceleration, unsigned MaxOrder>
class field_type_accelerator_pool :
	public pool<field_type_accelerator<Field, Family, acceleration::hard>, MaxOrder>
{
};


template <class Field, class Family, unsigned MaxOrder>
class field_type_accelerator_pool<Field, Family, acceleration::soft, MaxOrder> :
	public pool<typename quadrature_type<Family, typename Field::domain_t>::type, MaxOrder>
{
public:
	typedef pool<typename quadrature_type<Family, typename Field::domain_t>::type, MaxOrder> base_t;
	typedef field_type_accelerator<Field, Family, acceleration::soft> accelerator_t;

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

