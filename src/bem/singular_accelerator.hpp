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

/**
* \brief two iterators that can traverse two parallel containers
* \tparam Iter type of the iterator
* \todo should be replaced by matrix_iterator
*/
template <class Iter1, class Iter2 = Iter1>
class dual_iterator : private std::pair<Iter1, Iter2>
{
public:
	/** \brief the base type for abbreviations */
	typedef std::pair<Iter1, Iter2> base_t;

	/**
	* \brief constructor initialising all members
	* \param [in] prime the outer iterator
	* \param [in] sec the internal iterator
	*/
	dual_iterator(Iter1 prime, Iter2 sec)
		: base_t(prime, sec)
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
	typename Iter1::value_type const &get_prime(void) const
	{
		return *base_t::first;
	}

	/**
	* \brief return the secondary iterator's pointed value
	* \return the secondary iterator's pointed value
	*/
	typename Iter2::value_type const &get_sec(void) const
	{
		return *base_t::second;
	}
};


/**
* \brief a dual iterator to store a test and a trial quadrature iterator
* \tparam quadrature_iterator_t the iterator type of the quadratures
*/
template <class test_iterator_t, class trial_iterator_t>
class singular_quadrature_iterator : public dual_iterator<test_iterator_t, trial_iterator_t>
{
public:
	/** \brief the base type */
	typedef dual_iterator<test_iterator_t, trial_iterator_t> base_t;
	/** \brief the value type of both quadratures */
	typedef typename test_iterator_t::value_type value_type;

	/** \brief constructor initialising members
	* \param [in] test the test quadrature
	* \param [in] trial the trial quadrature
	*/
	singular_quadrature_iterator(test_iterator_t test, trial_iterator_t trial)
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
*/
template <class Formalism, class Kernel, class TestField, class TrialField>
class singular_accelerator;


/**
* \brief specialisation of ::singular_accelerator for the general formalism
* \tparam Kernel the kernel that is integrated
* \tparam TestField the test field over which integration is performed
* \tparam TrialField the trial field over which integration is performed
*/
template <class Kernel, class TestField, class TrialField>
class singular_accelerator<formalism::general, Kernel, TestField, TrialField>
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
	typedef singular_quadrature_iterator<
		typename test_quadrature_t::const_iterator,
		typename trial_quadrature_t::const_iterator> iterator;

	/** \brief indicates whether ::FACE_MATCH is possible */
	static const bool face_match_possible = std::is_same<test_field_t, trial_field_t>::value;
	/** \brief the singular quadrature order required by the kernel */
	static unsigned const singular_quadrature_order = kernel_traits<kernel_t>::singular_quadrature_order;

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
	iterator begin(void) const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return iterator(m_face_test_quadrature.begin(),
				m_face_trial_quadrature.begin());
			break;
		case EDGE_MATCH:
			return iterator(
				m_edge_test_quadrature[m_cur_overlap.get_ind1()].begin(),
				m_edge_trial_quadrature[m_cur_overlap.get_ind2()].begin());
			break;
		case CORNER_MATCH:
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()].begin(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()].begin());
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
	iterator end() const
	{
		switch (m_sing_type)
		{
		case FACE_MATCH:
			return iterator(
				m_face_test_quadrature.end(),
				m_face_trial_quadrature.end());
			break;
		case EDGE_MATCH:
			return iterator(
				m_edge_test_quadrature[m_cur_overlap.get_ind1()].end(),
				m_edge_trial_quadrature[m_cur_overlap.get_ind2()].end());
			break;
		case CORNER_MATCH:
			return iterator(
				m_corner_test_quadrature[m_cur_overlap.get_ind1()].end(),
				m_corner_trial_quadrature[m_cur_overlap.get_ind2()].end());
			break;
		case REGULAR:
			throw("Cannot return singular quadrature for regular type");
			break;
		default:
			throw("Unknown singularity type");
		}
	}

private:
	void generate_face(std::false_type)
	{
	}

	void generate_face(std::true_type)
	{
		// construct facials
		quad_factory_t::template generate<FACE_MATCH>(
			m_face_test_quadrature,
			m_face_trial_quadrature,
			singular_quadrature_order);
	}

public:
	/**
	* \brief constructor allocating space for the quadratures
	*/
	singular_accelerator(void)
	{
		// generate face quadratures with separate function
		// so that it does not compile of not needed (different domains)
		generate_face(std::integral_constant<bool, face_match_possible>());

		// create temporary quadratures for rotating
		test_quadrature_t test_edge_q, test_corner_q;
		trial_quadrature_t trial_edge_q, trial_corner_q;

		// create test quadrature singular on the first corner
		quad_factory_t::template generate<CORNER_MATCH>(
			test_corner_q, trial_corner_q, singular_quadrature_order);
		// create test quadrature singular on the first edge
		quad_factory_t::template generate<EDGE_MATCH>(
			test_edge_q, trial_edge_q, singular_quadrature_order);

		// rotate test quadratures
		for	(unsigned i = 0; i < n_test_corners; ++i)
		{
			// fill transform coordinates
			Eigen::Matrix<scalar_t, n_test_corners, test_domain_t::dimension> test_corners;
			for (unsigned k = 0; k < n_test_corners; ++k)
				test_corners.row(k) = test_domain_t::get_corner((i+k) % n_test_corners);

			// rotate
			m_corner_test_quadrature[i] =
				test_corner_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners);
			m_edge_test_quadrature[i] =
				test_edge_q.template transform<isoparam_shape_set<test_domain_t> >(test_corners);
		}
		// rotate trial quads
		for (unsigned i = 0; i < n_trial_corners; ++i)
		{
			// fill transform coordinates
			Eigen::Matrix<scalar_t, n_trial_corners, trial_domain_t::dimension> trial_corners;
			for (unsigned k = 0; k < n_trial_corners; ++k)
				trial_corners.row(k) = trial_domain_t::get_corner((i+k) % n_trial_corners);

			// rotate
			m_corner_trial_quadrature[i] =
				trial_corner_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners);
				
			// when dealing with the EDGE_MATCH case, we need to take into consideration that the test and
			// trial elements contain the singular edge in opposite nodal order. Therefore, if the
			// test element's singular edge is 0-1, then the trial element's singular edge should be 1-0
			// this numbering is implemented below: i+1-k instead of i+k

			// fill transform coordinates
			for (unsigned k = 0; k < n_trial_corners; ++k)
				trial_corners.row(k) =
				trial_domain_t::get_corner((i+1+n_trial_corners-k) % n_trial_corners);

			// rotate but keel weights unadjusted
			m_edge_trial_quadrature[i] = trial_edge_q.template transform<isoparam_shape_set<trial_domain_t> >(trial_corners);
			m_edge_trial_quadrature[i] *= -1.0;
		}
	}

protected:
	singularity_type m_sing_type;	/**< \brief the current singulartity type */
	element_overlapping m_cur_overlap;	/**< \brief the current overlapping state */

	/**\brief singular quadrature on the test elem for FACE_MATCH case */
	test_quadrature_t m_face_test_quadrature;
	/**\brief singular quadrature on the trial elem for FACE_MATCH case */
	trial_quadrature_t m_face_trial_quadrature;

	/**\brief singular quadratures on the test elem for CORNER_MATCH case */
	test_quadrature_t m_corner_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for CORNER_MATCH case */
	trial_quadrature_t m_corner_trial_quadrature[trial_domain_t::num_corners];

	/**\brief singular quadratures on the test elem for EDGE_MATCH case */
	test_quadrature_t m_edge_test_quadrature[test_domain_t::num_corners];
	/**\brief singular quadratures on the trial elem for EDGE_MATCH case */
	trial_quadrature_t m_edge_trial_quadrature[trial_domain_t::num_corners];
};



/**
* \brief specialisation the singular accelerator for the collocational case
* \tparam Kernel the kernel to integrate
* \tparam TestElem the test element type
* \tparam TrialField the trial field type
*/
template <class Kernel, class TestField, class TrialField>
class singular_accelerator<formalism::collocational, Kernel, TestField, TrialField>
{
	// CRTP check
	static_assert(std::is_base_of<kernel_base<Kernel>, Kernel>::value,
		"The kernel field must be derived from kernel_base<Kernel>");
	static_assert(std::is_base_of<field_base<TrialField>, TrialField>::value,
		"The trial field must be derived from field_base<TrialField>");
public:
	/** \brief template argument as nested type */
	typedef Kernel kernel_t;
	/** \brief the test field type */
	typedef TestField test_field_t;
	/** \brief template argument as nested type */
	typedef TrialField trial_field_t;

	typedef typename test_field_t::elem_t test_elem_t;	/**< \brief the test elem type */
	typedef typename trial_field_t::elem_t trial_elem_t;	/**< \brief trial elem type */

	typedef typename test_elem_t::domain_t test_domain_t;	/**< \brief test domain type */
	typedef typename trial_elem_t::domain_t trial_domain_t;	/**< \brief trial domain type */

	/** \brief trial L-set type */
	typedef typename trial_elem_t::lset_t trial_lset_t;
	/** \brief test N-set type */
	typedef typename test_field_t::nset_t test_nset_t;

	/** \brief quadrature family */
	typedef typename kernel_traits<kernel_t>::quadrature_family_t quadrature_family_t;
	/** \brief the trial quadrature type */
	typedef typename quadrature_type<	
		quadrature_family_t, trial_domain_t
	>::type trial_quadrature_t;
	/** \brief quadrature element type (it should be the same for test and trial) */
	typedef typename trial_quadrature_t::quadrature_elem_t quadrature_elem_t;
	/** \brief the singular quadrature order required by the kernel */
	static unsigned const singular_quadrature_order = kernel_traits<kernel_t>::singular_quadrature_order;

	/** \brief the duffy quadrature type */
	typedef duffy_quadrature<quadrature_family_t, trial_lset_t> trial_duffy_t;

	/** \brief iterator type of the singular quadrature */
	typedef typename trial_quadrature_t::const_iterator iterator;


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
	trial_quadrature_t const &get_trial_quadrature(unsigned idx) const
	{
		return m_face_trial_quadratures[idx];
	}

	/**
	* \brief constructor allocating memory for the quadratures
	*/
	singular_accelerator(void)
	{
		if (face_match_possible)
		{
			for (unsigned idx = 0; idx < test_nset_t::num_nodes; ++idx)
				m_face_trial_quadratures[idx] += trial_duffy_t::on_face(
					singular_quadrature_order, test_nset_t::corner_at(idx));
		}
	}

protected:
	singularity_type m_sing_type;	/**< \brief the actual singularity type */
	/** \brief pointer to the trial singular quadrature */
	trial_quadrature_t m_face_trial_quadratures[test_nset_t::num_nodes];	
};

#endif // SINGULAR_ACCELERATOR_HPP_INCLUDED

