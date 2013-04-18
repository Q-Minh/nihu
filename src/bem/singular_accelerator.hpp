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

	dual_iterator(Primary prime, Secondary sec, unsigned Nsec, dual_iterator_mode mode)
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

	typename Primary::value_type const &get_prime(void) const
	{
		return *m_prime;
	}

	typename Secondary::value_type const &get_sec(void) const
	{
		return *m_sec;
	}

protected:
	Primary m_prime;
	Secondary m_sec;
	unsigned const m_Nsec;
	unsigned m_cntr;
	dual_iterator_mode m_mode;
	Secondary m_sec_start;
};

template <class quadrature_iterator_t>
class singular_quadrature_iterator : public dual_iterator<quadrature_iterator_t, quadrature_iterator_t>
{
public:
	typedef dual_iterator<quadrature_iterator_t, quadrature_iterator_t> base_t;

	singular_quadrature_iterator(quadrature_iterator_t test, quadrature_iterator_t trial, unsigned Ntrial, typename base_t::dual_iterator_mode mode)
		: base_t(test, trial, Ntrial, mode)
	{
	}

	typename quadrature_iterator_t::value_type const &get_test_quadrature_elem(void) const
	{
		return this->get_prime();
	}

	typename quadrature_iterator_t::value_type const &get_trial_quadrature_elem(void) const
	{
		return this->get_sec();
	}
};

enum singularity_type {
	REGULAR,
	FACE_MATCH,
	EDGE_MATCH,
	CORNER_MATCH
};

/**
 * \brief an accelerator class that stores singular quadratures for different singularity types
 * \tparam Kernel the kernel that is integrated
 * \tparam TestField the test field over which integration is performed
 * \tparam TrialField the trial field over which integration is performed
 */
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
	typedef Kernel kernel_t;	/**< \brief template argument as nested type */
	typedef TestField test_field_t;	/**< \brief template argument as nested type */
	typedef TrialField trial_field_t;	/**< \brief template argument as nested type */

	typedef typename test_field_t::elem_t test_elem_t;	/**< \brief the test elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;	/**< \brief trial elem type */

	typedef typename quadrature_type<
		typename kernel_traits<kernel_t>::quadrature_family_t,
		typename test_elem_t::domain_t
	>::type test_quadrature_t;			/**< \brief the test quadrature type */
	typedef typename quadrature_type<
		typename kernel_traits<kernel_t>::quadrature_family_t,
		typename trial_elem_t::domain_t
	>::type trial_quadrature_t;			/**< \brief the trial quadrature type */

	/**\brief the quadrature element type (it should be the same for test and trial) */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;

	/**\brief the dual iterator type of te singular quadrature */
	typedef singular_quadrature_iterator<typename test_quadrature_t::iterator> iterator;

	/**\brief indicates whether FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_elem_t, trial_elem_t>::value;

	singularity_type eval(test_field_t const &test_field, trial_field_t const &trial_field) const
	{
		// check face match for same element types
		if (face_match_possible)
			if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
				return FACE_MATCH;

		/*
		element_overlapping eo = test_field.get_elem().get_overlapping(trial_field.get_elem());
		// check edge match
		if (eo.get_num() > 1)
			return EDGE_MATCH;
		// check corner match
		if (eo.get_num() > 0)
			return CORNER_MATCH;
		*/
		return REGULAR;
	}

	iterator cbegin(singularity_type sing_type) const
	{
		switch (sing_type)
		{
		case FACE_MATCH:
			return iterator(
				m_face_test_quadrature->begin(),
				m_face_trial_quadrature->begin(),
				m_face_trial_quadrature->size() / m_face_test_quadrature->size(),
				iterator::base_t::CONTINUOUS);
			break;
		case EDGE_MATCH:
			throw("Unimplemented singular quadrature");
			break;
		case CORNER_MATCH:
			throw("Unimplemented singular quadrature");
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

	iterator cend(singularity_type sing_type) const
	{
		switch (sing_type)
		{
		case FACE_MATCH:
			return iterator(
				m_face_test_quadrature->end(),
				m_face_trial_quadrature->end(),
				m_face_trial_quadrature->size() / m_face_test_quadrature->size(),
				iterator::CONTINUOUS);
			break;
		case EDGE_MATCH:
			throw("Unimplemented singular quadrature");
			break;
		case CORNER_MATCH:
			throw("Unimplemented singular quadrature");
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

	singular_accelerator(void)
	{
		if (face_match_possible)
		{
			m_face_test_quadrature = new test_quadrature_t(9);
			m_face_trial_quadrature = new trial_quadrature_t();
			for (auto it = m_face_test_quadrature->begin();
				it != m_face_test_quadrature->end(); ++it)
				*m_face_trial_quadrature +=
				trial_quadrature_t::singular_quadrature_inside(9, it->get_xi());
		}
		else
		{
			m_face_test_quadrature = NULL;
			m_face_trial_quadrature = NULL;
		}
	}

	~singular_accelerator()
	{
		if (face_match_possible)
		{
			delete m_face_test_quadrature;
			delete m_face_trial_quadrature;
		}
	}

protected:
	/** \todo The code does not compile if the quadratures are stored statically,
	 * not allocated dynamically. And we do not understand this sickness.
	 */
	/**\brief singular quadrature on the test elem for FACE_MATCH case */
	test_quadrature_t *m_face_test_quadrature;
	/**\brief singular quadrature on the trial elem for FACE_MATCH case */
	trial_quadrature_t *m_face_trial_quadrature;
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
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;
	typedef typename quadrature_type<
		typename kernel_traits<kernel_t>::quadrature_family_t,
		trial_domain_t
	>::type trial_quadrature_t;

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

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

