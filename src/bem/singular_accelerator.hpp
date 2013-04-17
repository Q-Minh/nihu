/**
 * \file singular_accelerator.hpp
 * \brief implementation of class ::singular_accelerator
 * \author Peter Fiala and Peter Rucz
 */

 #ifndef SINGULAR_ACCELERATOR_HPP_INCLUDED
#define SINGULAR_ACCELERATOR_HPP_INCLUDED

#include "field.hpp"


template <class Primary, class Secondary>
class dual_iterator
{
public:
	enum dual_iterator_mode { CIRCULAR, CONTINUOUS };

	dual_iterator(Primary const *prime, Secondar const *sec, unsigned Nsec, dual_iterator_mode mode)
		: m_prime(prime), m_sec(sec), m_Nsec(Nsec), m_cntr(0), m_mode(mode), m_sec_start(m_sec)
	{
	}

	dual_iterator &operator++(void)
	{
		++m_sec;
		++m_cntr;
		if (m_cntr == m_Nsec)
		{
			m_cntr = 0;
			++m_prime;
			if (m_mode == CIRCULAR)
				m_sec = m_sec_start;
		}
		return *this;
	}

	bool operator!=(dual_iterator const &other)
	{
		return m_prime != other.m_prime || m_sec != other.m_sec;
	}

	Primary const &get_prime(void) const
	{
		return *m_prime;
	}

	Secondary const &get_sec(void) const
	{
		return *m_sec;
	}

protected:
	dual_iterator_mode m_mode;
	Primary const *m_prime;
	Secondary const *m_sec;
	Secondary const *m_sec_start;
	unsigned const m_Nsec;
	unsigned m_cntr;
};


enum singularity_type {
	REGULAR,
	FACE_MATCH,
	EDGE_MATCH,
	CORNER_MATCH
};

/// The general case
template <class Kernel, class TestField, class TrialField>
class singular_accelerator
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"The kernel field must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TestField>, TestField>::value,
		"The test field must be derived from field_base<TestField>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	typedef Kernel kernel_t;
	typedef TestField test_field_t;
	typedef TrialField trial_field_t;
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	typedef typename test_field_t::elem_t::domain_t test_domain_t;
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;

	class iterator : public dual_iterator<quadrature_elem_t, quadrature_elem_t>
	{
		typedef dual_iterator<quadrature_elem_t, quadrature_elem_t> base_t;
	public:
		iterator(quadrature_elem_t const *test, quadrature_elem_t const *trial, unsigned Ntrial, dual_iterator::dual_iterator_mode mode)
			: base_t(test, trial, Ntrial, mode)
		{
		}

		quadrature_elem_t const &get_test_quadrature_elem(void) const
		{
			return this->get_prime();
		}

		quadrature_elem_t const &get_trial_quadrature_elem(void) const
		{
			return this->get_sec();
		}
	};

	singularity_type eval(TestField const &, TrialField const &) const
	{
		return REGULAR;
	}

	test_quadrature_t const *get_test_quadrature(void) const
	{
		return NULL;
	}

	trial_quadrature_t const *get_trial_quadrature(void) const
	{
		return NULL;
	}
};


/// The collocational constant case
template <class Kernel, class Elem, class TrialFieldOption, class TrialDiracOption>
class singular_accelerator<
	Kernel,
	field<Elem, constant_field, dirac_field>,
	field<Elem, TrialFieldOption, TrialDiracOption> >
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");
public:
	typedef Kernel kernel_t;
	typedef field<Elem, constant_field, dirac_field> test_field_t;
	typedef field<Elem, TrialFieldOption, TrialDiracOption> trial_field_t;
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;

	singular_accelerator()
		: m_quadrature(trial_quadrature_t::singular_quadrature_inside(9, trial_domain_t::get_center()))
	{
	}

	singularity_type eval(test_field_t const &test_field, trial_field_t const &trial_field) const
	{
		if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			return FACE_MATCH;
		return REGULAR;
	}

	trial_quadrature_t const *get_trial_quadrature(void) const
	{
		return &m_quadrature;
	}

protected:
	trial_quadrature_t const m_quadrature;
};


/// The general galerkin case
template <class Kernel, class Elem, class TestOption, class TrialField>
class singular_accelerator<
	Kernel,
	field<Elem, TestOption, function_field>,
	TrialField>
{
	// CRTP test
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"Kernel must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"TrialField must be derived from field_base<TrialField>");
public:
	typedef field<Elem, TestOption, function_field> test_field_t;
	typedef TrialField  trial_field_t;
	typedef typename test_field_t::elem_t test_elem_t;
	typedef typename trial_field_t::elem_t trial_elem_t;

	singularity_type eval(test_field_t const &test_field, trial_field_t const &trial_field) const
	{
		// check face match for same element types
		if (std::is_same<test_elem_t, trial_elem_t>::value)
			if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
				return FACE_MATCH;

		element_overlapping eo = test_field.get_elem().get_overlapping(trial_field.get_elem());
		// check edge match
		if (eo.get_num() > 1)
			return EDGE_MATCH;
		// check corner match
		if (eo.get_num() > 0)
			return CORNER_MATCH;
		return REGULAR;
	}
};

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

