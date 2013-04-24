/**
* \file singular_accelerator.hpp
* \brief implementation of class ::singular_accelerator
* \author Peter Fiala and Peter Rucz
*/

#ifndef SINGULAR_ACCELERATOR_HPP_INCLUDED
#define SINGULAR_ACCELERATOR_HPP_INCLUDED

#include "field.hpp"
#include "duffy_quadrature.hpp"

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

	typedef typename test_elem_t::domain_t test_domain_t;
	typedef typename trial_elem_t::domain_t trial_domain_t;

	typedef typename test_elem_t::lset_t test_lset_t;
	typedef typename trial_elem_t::lset_t trial_lset_t;
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;

	/** \brief the test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;

	/**\brief the quadrature element type (it should be the same for test and trial) */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;

	typedef duffy_quadrature<quadrature_family_t, test_lset_t> test_duffy_t;
	typedef duffy_quadrature<quadrature_family_t, trial_lset_t> trial_duffy_t;

	/**\brief the dual iterator type of te singular quadrature */
	typedef singular_quadrature_iterator<typename test_quadrature_t::iterator> iterator;

	/**\brief indicates whether FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_elem_t, trial_elem_t>::value;

	bool is_singular(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		// check face match for same element types
		if (face_match_possible)
		{
			if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			{
				m_sing_type = FACE_MATCH;
				return true;
			}
		}

		m_cur_overlap = test_field.get_elem().get_overlapping(trial_field.get_elem());
//		// check edge match
//		if (eo.get_num() > 1)
//			return EDGE_MATCH;
		// check corner match
		if (m_cur_overlap.get_num() == 1)
		{
			m_sing_type = CORNER_MATCH;
			return true;
		}

		m_sing_type = REGULAR;
		return false;
	}

	iterator cbegin(void) const
	{
		switch (m_sing_type)
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
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()]->begin(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->begin(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->size(),
				iterator::base_t::CIRCULAR);
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

	iterator cend() const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return iterator(
				m_face_test_quadrature->end(),
				m_face_trial_quadrature->end(),
				m_face_trial_quadrature->size() / m_face_test_quadrature->size(),
				iterator::base_t::CONTINUOUS);
			break;
		case EDGE_MATCH:
			throw("Unimplemented singular quadrature");
			break;
		case CORNER_MATCH:
			/** \todo It is very dangerous that we need to handle the CIRCULAR case with { end-begin } dual quadrature */
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()]->end(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->begin(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->size(),
				iterator::base_t::CIRCULAR);
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
		/** \todo kernel should tell the singularity order */
		unsigned const SINGULARITY_ORDER = 9;

		if (face_match_possible)
		{
			m_face_test_quadrature = new test_quadrature_t(SINGULARITY_ORDER);
			m_face_trial_quadrature = new trial_quadrature_t();
			for (auto it = m_face_test_quadrature->begin();
				it != m_face_test_quadrature->end(); ++it)
				*m_face_trial_quadrature += trial_duffy_t::on_face(SINGULARITY_ORDER, it->get_xi());
		}
		else
		{
			m_face_test_quadrature = NULL;
			m_face_trial_quadrature = NULL;
		}

		for (unsigned i = 0; i < test_domain_t::num_corners; ++i)
		{
			m_corner_test_quadrature[i] = new test_quadrature_t();
			*m_corner_test_quadrature[i] += test_duffy_t::on_corner(SINGULARITY_ORDER, i);
		}
		for (unsigned i = 0; i < trial_domain_t::num_corners; ++i)
		{
			m_corner_trial_quadrature[i] = new trial_quadrature_t();
			*m_corner_trial_quadrature[i] += trial_duffy_t::on_corner(SINGULARITY_ORDER, i);
		}
	}

	~singular_accelerator()
	{
		if (face_match_possible)
		{
			delete m_face_test_quadrature;
			delete m_face_trial_quadrature;
		}

		for (unsigned i = 0; i < test_domain_t::num_corners; ++i)
			delete m_corner_test_quadrature[i];
		for (unsigned i = 0; i < trial_domain_t::num_corners; ++i)
			delete m_corner_trial_quadrature[i];
	}

protected:
	singularity_type m_sing_type;
	element_overlapping m_cur_overlap;

	/** \todo The code does not compile if the quadratures are stored statically,
	 * not allocated dynamically. And we do not understand this sickness.
	 */

	/**\brief singular quadrature on the test elem for FACE_MATCH case */
	test_quadrature_t *m_face_test_quadrature;
	/**\brief singular quadrature on the trial elem for FACE_MATCH case */
	trial_quadrature_t *m_face_trial_quadrature;

	test_quadrature_t *m_corner_test_quadrature[test_domain_t::num_corners];
	trial_quadrature_t *m_corner_trial_quadrature[trial_domain_t::num_corners];
};



template <class Kernel, class TestElem, class TrialField>
class singular_accelerator<Kernel, field<TestElem, constant_field, dirac_field>, TrialField>
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"The kernel field must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	typedef Kernel kernel_t;	/**< \brief template argument as nested type */
	typedef field<TestElem, constant_field, dirac_field> test_field_t;	/**< \brief template argument as nested type */
	typedef TrialField trial_field_t;	/**< \brief template argument as nested type */

	typedef typename test_field_t::elem_t test_elem_t;	/**< \brief the test elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;	/**< \brief trial elem type */

	typedef typename test_elem_t::domain_t test_domain_t;
	typedef typename trial_elem_t::domain_t trial_domain_t;

	typedef typename test_elem_t::lset_t test_lset_t;
	typedef typename trial_elem_t::lset_t trial_lset_t;
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;

	/** \brief the test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;

	/**\brief the quadrature element type (it should be the same for test and trial) */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;

	typedef duffy_quadrature<quadrature_family_t, test_lset_t> test_duffy_t;
	typedef duffy_quadrature<quadrature_family_t, trial_lset_t> trial_duffy_t;

	/**\brief the dual iterator type of te singular quadrature */
	typedef singular_quadrature_iterator<typename test_quadrature_t::iterator> iterator;

	/**\brief indicates whether FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_elem_t, trial_elem_t>::value;

	bool is_singular(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		if (face_match_possible)
		{
			if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			{
				m_sing_type = FACE_MATCH;
				return true;
			}
		}
		return false;
	}

	typename trial_quadrature_t::iterator cbegin(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return m_face_trial_quadrature->begin();
		}
	}

	typename trial_quadrature_t::iterator cend(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return m_face_trial_quadrature->end();
		}
	}

	singular_accelerator(void)
	{
		/** \todo kernel should tell the singularity order */
		unsigned const SINGULARITY_ORDER = 9;

		if (face_match_possible)
		{
			m_face_trial_quadrature = new trial_quadrature_t();
			*m_face_trial_quadrature += trial_duffy_t::on_face(SINGULARITY_ORDER, trial_domain_t::get_center());
		}
		else
		{
			m_face_trial_quadrature = NULL;
		}
	}

	~singular_accelerator()
	{
		if (face_match_possible)
		{
			delete m_face_trial_quadrature;
		}
	}

protected:
	singularity_type m_sing_type;
	trial_quadrature_t *m_face_trial_quadrature;
};



#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

