/**
* \file singular_accelerator.hpp
* \brief implementation of class ::singular_accelerator
* \author Peter Fiala and Peter Rucz
*/

#ifndef SINGULAR_ACCELERATOR_HPP_INCLUDED
#define SINGULAR_ACCELERATOR_HPP_INCLUDED

#include "singular_galerkin_quadrature.hpp"	// for the galerkin case
#include "duffy_quadrature.hpp"		// for the collocational case
#include "field.hpp"
#include <utility> // pair

/**
* \brief two iterators that can traverse two parallel containers
* \tparam Iter type of the iterator
*/
template <class Iter>
class dual_iterator : private std::pair<Iter, Iter>
{
public:
	/** \brief the base type for abbreviations */
	typedef std::pair<Iter, Iter> base_t;
	
	/**
	* \brief constructor initialising all members
	* \param [in] prime the outer iterator
	* \param [in] sec the internal iterator
	*/
	dual_iterator(Iter prime, Iter sec)
		: std::pair<Iter, Iter>(prime, sec)
	{
	}

	/**
	* \brief preincrement operator
	* \return reference to the incremented iterator
	*/
	dual_iterator &operator++(void)
	{
		++base_t::first;
		++base_t::second;
		return *this;
	}

	/**
	* \brief not equal operator
	* \param [in] other the other iterator
	* \return true if the iterators are different
	*/
	bool operator!=(dual_iterator const &other)
	{
		return base_t::first != other.base_t::first || base_t::second != other.base_t::second;
	}

	/**
	* \brief return the primary iterator's pointed value
	* \return the primary iterator's pointed value
	*/
	typename Iter::value_type const &get_prime(void) const
	{
		return *base_t::first;
	}

	/**
	* \brief return the secondary iterator's pointed value
	* \return the secondary iterator's pointed value
	*/
	typename Iter::value_type const &get_sec(void) const
	{
		return *base_t::second;
	}
};


/**
* \brief a dual iterator to store a test and a trial quadrature iterator
* \tparam quadrature_iterator_t the iterator type of the quadratures
*/
template <class quadrature_iterator_t>
class singular_quadrature_iterator : public dual_iterator<quadrature_iterator_t>
{
public:
	/** \brief the base type */
	typedef dual_iterator<quadrature_iterator_t> base_t;
	/** \brief the value type of both quadratures */
	typedef typename quadrature_iterator_t::value_type value_type;

	/** \brief constructor initialising members
	 * \param [in] test the test quadrature
	 * \param [in] trial the trial quadrature
	 */
	singular_quadrature_iterator(
		quadrature_iterator_t test,
		quadrature_iterator_t trial)
		: base_t(test, trial)
	{
	}

	/** \brief return the test quadrature element
	* \return test quadrature element
	*/
	value_type const &get_test_quadrature_elem(void) const
	{
		return this->get_prime();
	}

	/** \brief return the trial quadrature element
	* \return trial quadrature element
	*/
	value_type const &get_trial_quadrature_elem(void) const
	{
		return this->get_sec();
	}
};

/**
* \brief an accelerator class that stores singular quadratures for different singularity types
* \tparam Kernel the kernel that is integrated
* \tparam TestField the test field over which integration is performed
* \tparam TrialField the trial field over which integration is performed
* \todo templating this class on the kernel is sick.
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

	typedef typename test_field_t::elem_t::domain_t test_domain_t;		/**< \brief test domain */
	typedef typename trial_field_t::elem_t::domain_t trial_domain_t;	/**< \brief trial domain */

	/** \brief number of test domain corners */
	static unsigned const n_test_corners = test_domain_t::num_corners;
	/** \brief number of trial domain corners */
	static unsigned const n_trial_corners = trial_domain_t::num_corners;

	/** \brief the scalar type of the quadrature */
	typedef typename test_domain_t::scalar_t scalar_t;
	// report compiler error if the trial and test scalar types are different
	static_assert(std::is_same<scalar_t, typename trial_domain_t::scalar_t>::value,
		"The test and trial domain scalar types must match");

	/** \brief the quadrature family obtained from the kernel */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;

	/** \brief the singular galerkin quadrature generator class */
	typedef singular_galerkin_quadrature<quadrature_family_t, test_domain_t, trial_domain_t> quad_factory_t;

	/** \brief the test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;

	/** \brief the quadrature element type (it should be the same for test and trial) */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;
	// report compiler error if the trial and test quadrature elements are different
	static_assert(std::is_same<quadrature_elem_t, typename trial_quadrature_t::quadrature_elem_t>::value,
		"The trial and test quadrature elements must be of the same type");

	/** \brief the dual iterator type of the singular quadrature */
	typedef singular_quadrature_iterator<typename test_quadrature_t::iterator> iterator;

	/** \brief indicates whether ::FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_field_t, trial_field_t>::value;

	/**
	* \brief determine if singular integration is needed and store singularity type
	* \param [in] test_field the test field
	* \param [in] trial_field the trial_field
	* \return true if singular integration is needed
	*/
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
		// check edge match
		if (m_cur_overlap.get_num() > 1)
		{
			m_sing_type = EDGE_MATCH;
			return true;
		}
		// check corner match
		if (m_cur_overlap.get_num() == 1)
		{
			m_sing_type = CORNER_MATCH;
			return true;
		}

		m_sing_type = REGULAR;
		return false;
	}

	/**
	* \brief return begin iterator of the singular quadrature
	* \return begin iterator of the singular quadrature
	*/
	iterator cbegin(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return iterator(m_face_test_quadrature->begin(),
				m_face_trial_quadrature->begin());
			break;
		case EDGE_MATCH:
			return iterator(
				m_edge_test_quadrature[m_cur_overlap.get_ind1()]->begin(),
				m_edge_trial_quadrature[m_cur_overlap.get_ind2()]->begin());
			break;
		case CORNER_MATCH:
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()]->begin(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->begin());
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

	/**
	* \brief return end iterator of the singular quadrature
	* \return end iterator of the singular quadrature
	*/
	iterator cend() const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return iterator(
				m_face_test_quadrature->end(),
				m_face_trial_quadrature->end());
			break;
		case EDGE_MATCH:
			return iterator(
				m_edge_test_quadrature[m_cur_overlap.get_ind1()]->end(),
				m_edge_trial_quadrature[m_cur_overlap.get_ind2()]->end());
			break;
		case CORNER_MATCH:
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()]->end(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()]->end());
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

	/**
	* \brief constructor allocating space for the quadratures
	*/
	singular_accelerator(void)
		: m_face_test_quadrature(NULL), m_face_trial_quadrature(NULL)
	{
		/** \todo kernel should tell the singularity order */
		unsigned const SINGULARITY_ORDER = 2;

		if (face_match_possible)
		{
			// construct facials
			m_face_test_quadrature = new test_quadrature_t;
			m_face_trial_quadrature = new trial_quadrature_t;
			quad_factory_t::template generate<FACE_MATCH>(*m_face_test_quadrature,
				*m_face_trial_quadrature, SINGULARITY_ORDER);
		}

		Eigen::Matrix<scalar_t, n_test_corners, test_domain_t::dimension> test_corners;
		Eigen::Matrix<scalar_t, n_trial_corners, trial_domain_t::dimension> trial_corners;

		// create temporary quadratures for rotating
		test_quadrature_t test_edge_q, test_corner_q;
		trial_quadrature_t trial_edge_q, trial_corner_q;

		// create quads singular on the first edge
		quad_factory_t::template generate<EDGE_MATCH>(
			test_edge_q, trial_edge_q, SINGULARITY_ORDER);
		// create quads singular on the first corner
		quad_factory_t::template generate<CORNER_MATCH>(
			test_corner_q, trial_corner_q, SINGULARITY_ORDER);

		// rotate test quads
		for	(unsigned i = 0; i < n_test_corners; ++i)
		{
			// fill transform coordinates
			for (unsigned k = 0; k < n_test_corners; ++k)
				test_corners.row(k) = test_domain_t::get_corners()[(i+k) % n_test_corners].transpose();

			// rotate
			m_corner_test_quadrature[i] = new test_quadrature_t(
				test_corner_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners));
			m_edge_test_quadrature[i] = new test_quadrature_t(
				test_edge_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners));
		}
		// rotate trial quads
		for (unsigned i = 0; i < n_trial_corners; ++i)
		{
			// fill transform coordinates
			for (unsigned k = 0; k < n_trial_corners; ++k)
				trial_corners.row(k) = test_domain_t::get_corners()[(i+k) % n_trial_corners].transpose();

			// rotate
			m_corner_trial_quadrature[i] = new trial_quadrature_t(
				trial_corner_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners));
			m_edge_trial_quadrature[i] = new trial_quadrature_t(
				trial_edge_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners));
		}
	}

	/**
	* \brief destructor freeing quadratures
	*/
	~singular_accelerator()
	{
		if (face_match_possible)
		{
			delete m_face_test_quadrature;
			delete m_face_trial_quadrature;
		}

		for (unsigned i = 0; i < n_test_corners; ++i)
			delete m_corner_test_quadrature[i];
		for (unsigned i = 0; i < n_trial_corners; ++i)
			delete m_corner_trial_quadrature[i];

		for (unsigned i = 0; i < n_test_corners; ++i)
			delete m_edge_test_quadrature[i];
		for (unsigned i = 0; i < n_trial_corners; ++i)
			delete m_edge_trial_quadrature[i];
	}

protected:
	singularity_type m_sing_type;	/**< \brief the current singulartity type */
	element_overlapping m_cur_overlap;	/**< \brief the current overlapping state */

	/** \todo The code does not compile if the quadratures are stored statically,
	* not allocated dynamically. And we do not understand this sickness.
	*/

	/**\brief singular quadrature on the test elem for FACE_MATCH case */
	test_quadrature_t *m_face_test_quadrature;
	/**\brief singular quadrature on the trial elem for FACE_MATCH case */
	trial_quadrature_t *m_face_trial_quadrature;

	/**\brief singular quadratures on the test elem for CORNER_MATCH case */
	test_quadrature_t *m_corner_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for CORNER_MATCH case */
	trial_quadrature_t *m_corner_trial_quadrature[trial_domain_t::num_corners];

	/**\brief singular quadratures on the test elem for EDGE_MATCH case */
	test_quadrature_t *m_edge_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for EDGE_MATCH case */
	trial_quadrature_t *m_edge_trial_quadrature[trial_domain_t::num_corners];
};



/**
* \brief specialisation the singular accelerator for the collocational case
* \tparam Kernel the kernel to integrate
* \tparam TestElem the test element type
* \tparam TrialField the trial field type
*/
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
	/** \brief the test field type */
	typedef field<TestElem, constant_field, dirac_field> test_field_t;	
	typedef TrialField trial_field_t;	/**< \brief template argument as nested type */

	typedef typename test_field_t::elem_t test_elem_t;	/**< \brief the test elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;	/**< \brief trial elem type */

	typedef typename test_elem_t::domain_t test_domain_t;	/**< \brief test domain type */
	typedef typename trial_elem_t::domain_t trial_domain_t;	/**< \brief trial domain type */

	typedef typename trial_elem_t::lset_t trial_lset_t;	/**< \brief trial L-set type */

	/** \brief quadrature family */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<	
		quadrature_family_t, trial_domain_t
	>::type trial_quadrature_t;
	/** \brief quadrature element type (it should be the same for test and trial) */
	typedef typename trial_quadrature_t::quadrature_elem_t quadrature_elem_t;

	/** \brief the duffy quadrature type */
	typedef duffy_quadrature<quadrature_family_t, trial_lset_t> trial_duffy_t;

	/** \brief iterator type of the singular quadrature */
	typedef typename trial_quadrature_t::iterator iterator;


	/** \brief indicates whether FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_elem_t, trial_elem_t>::value;

	/**
	* \brief determine if singular integration is needed and store singularity type
	* \param [in] test_field the test field
	* \param [in] trial_field the trial_field
	* \return true if singular integration is needed
	*/
	bool is_singular(test_field_t const &test_field, trial_field_t const &trial_field)
	{
		if (face_match_possible && test_field.get_elem().get_id() == trial_field.get_elem().get_id())
		{
			m_sing_type = FACE_MATCH;
			return true;
		}
		return false;
	}

	/**
	* \brief return begin iterator of the singular quadrature
	* \return begin iterator of the singular quadrature
	*/
	iterator cbegin(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return m_face_trial_quadrature->begin();
		default:
			throw "unhandled quadrature type";
		}
	}

	/**
	* \brief return end iterator of the singular quadrature
	* \return end iterator of the singular quadrature
	*/
	iterator cend(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return m_face_trial_quadrature->end();
		default:
			throw "unhandled quadrature type";
		}
	}

	/**
	* \brief constructor allocating memory for the quadratures
	*/
	singular_accelerator(void)
	{
		/** \todo kernel should tell the singularity order */
		unsigned const SINGULARITY_ORDER = 9;

		if (face_match_possible)
		{
			m_face_trial_quadrature = new trial_quadrature_t();
			*m_face_trial_quadrature += trial_duffy_t::on_face(
				SINGULARITY_ORDER, trial_domain_t::get_center());
		}
		else
			m_face_trial_quadrature = NULL;
	}

	/**
	* \brief destructor freeing quadratures
	*/
	~singular_accelerator()
	{
		if (face_match_possible)
			delete m_face_trial_quadrature;
	}

protected:
	singularity_type m_sing_type;	/**< \brief the actual singularity type */
	/** \brief pointer to the trial singular quadrature */
	trial_quadrature_t *m_face_trial_quadrature;	
};

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

