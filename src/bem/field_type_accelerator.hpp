/**
 * \file field_type_accelerator.hpp
 * \ingroup quadrature
 * \brief declaration of class field_type_accelerator
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

#include "field.hpp"
#include "quadrature.hpp"
#include "../util/casted_iterator.hpp"
#include "../util/dual_range.hpp"


/** \brief collection of acceleration-types */
namespace acceleration
{
	/** \brief real acceleration */
	struct hard {};
	/** \brief view-acceleration */
	struct soft{};
}

template<class Field, class Family, class Acceleration>
class field_type_accelerator_elem;

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

	typename Field::nset_t::shape_t get_N(void) const
	{
		return Field::nset_t::eval_shape(base_t::get_xi());
	}
};

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

	field_type_accelerator_elem(base_t const &base) :
		base_t(base),
		m_nset(Field::nset_t::eval_shape(base_t::get_xi()))
	{
	}

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
	public quadrature_type<
		Family, typename Field::domain_t
	>::type
{
public:
	typedef typename quadrature_type<
		Family, typename Field::domain_t
	>::type base_t;

	typedef field_type_accelerator_elem<
		Field, Family, acceleration::soft
	> accelerator_elem_t;

	typedef casted_iterator<
		typename base_t::const_iterator,
		accelerator_elem_t
	> iterator;

	iterator begin(void) const
	{
		return base_t::begin();
	}

	iterator end(void) const
	{
		return base_t::end();
	}
};


struct index_t
{
	index_t(unsigned idx) :	m_idx(idx) {}
	unsigned m_idx;
};


template <class NSet>
class dirac_field_type_accelerator_elem :
	public index_t
{
public:
	constexpr typename NSet::scalar_t get_w(void) const
	{
		return 1.0;
	}

	constexpr typename NSet::shape_t get_N(void) const
	{
		return NSet::shape_t::Unit(index_t::m_idx);
	}

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

	bool operator !=(dirac_field_type_accelerator_iterator const &other)
	{
		return m_idx.m_idx != other.m_idx.m_idx;
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


template <class Field, class Family, class Acceleration>
class field_type_accelerator<Field, Family, Acceleration,
	typename std::enable_if<field_traits<Field>::is_dirac>::type
	>
{
public:
	typedef dirac_field_type_accelerator_elem<
		typename Field::nset_t
	> accelerator_elem_t;

	field_type_accelerator(void)
	{
	}

	typedef casted_iterator<
		dirac_field_type_accelerator_iterator,
		accelerator_elem_t
	> iterator;

	constexpr iterator begin(void) const
	{
		return dirac_field_type_accelerator_iterator(index_t(0));
	}

	constexpr iterator end(void) const
	{
		return dirac_field_type_accelerator_iterator(index_t(Field::num_nodes));
	}
};



template <class TestField, class TrialField, class IterationMode, class Family, class Acceleration>
class dual_field_type_accelerator :
	public dual_range<
		IterationMode,
		typename field_type_accelerator<TestField, Family, Acceleration>::iterator,
		typename field_type_accelerator<TrialField, Family, Acceleration>::iterator
	>
{
public:
	typedef typename field_type_accelerator<TestField, Family, Acceleration>::iterator test_iter_t;
	typedef typename field_type_accelerator<TrialField, Family, Acceleration>::iterator trial_iter_t;
	typedef dual_range<IterationMode, test_iter_t, trial_iter_t> base_t;

	dual_field_type_accelerator(
		test_iter_t const &test_accelerator,
		trial_iter_t const &trial_accelerator
	) : base_t(
		test_accelerator.begin(), test_accelerator.end(),
		trial_accelerator.begin(), trial_accelerator.end()
	)
	{
	}
};


/**
 * \brief container class to store accelerators with different quadrature orders
 * \tparam Field the field type that is accelerated
 */
template <class Field, class Family, class Acceleration>
class field_type_accelerator_pool :
	public std::vector<field_type_accelerator<Field, Family, Acceleration> *>
{
	// CRTP check
	static_assert(std::is_base_of<field_base<Field>, Field>::value,
		"Field must be derived from field_base<Field>");
public:
	/** \brief template argument as nested type */
	typedef Field field_t;
	/** \brief template argument as nested type */
	typedef Family family_t;

	/** \brief maximum order of quadratures */
	static const unsigned MAX_ORDER = 9;

	/** \brief allocate memory and initialise accelerators */
	field_type_accelerator_pool(void)
	{
		this->reserve(MAX_ORDER+1);
		for (unsigned order = 0; order <= MAX_ORDER; ++order)
			this->push_back(new field_type_accelerator<field_t, family_t, acceleration::hard>(order));
	}

	/** \brief free memory allocated for the accelerators */
	~field_type_accelerator_pool(void)
	{
		for (unsigned order = 0; order <= MAX_ORDER; ++order)
			delete (*this)[order];
	}
};

#endif // FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

