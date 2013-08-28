#ifndef ELEMENT_MATCH_HPP_INCLUDED
#define ELEMENT_MATCH_HPP_INCLUDED

#include "element.hpp"

/** \brief singularity types */
enum singularity_type {
	REGULAR,		/**< \brief no singularity */
	FACE_MATCH,	/**< \brief two elements are identical */
	EDGE_MATCH,	/**< \brief two elements share common edge */
	CORNER_MATCH	/**< \brief two elements share common corner */
};

namespace singularity
{
	typedef std::integral_constant<unsigned, REGULAR> regular_type;
	typedef std::integral_constant<unsigned, REGULAR> face_match_type;
	typedef std::integral_constant<unsigned, REGULAR> edge_match_type;
	typedef std::integral_constant<unsigned, REGULAR> corner_match_type;
}

/** \brief class describing the adjacency (match) state of two elements */
class element_match
{
public:
	/** \brief constructor
	 * \param [in] sing_type the singularity type
	 * \param [in] overlap the overlapping state
	 */
	element_match(
		singularity_type const &sing_type,
		element_overlapping const &overlap = element_overlapping()) :
		m_sing_type(sing_type), m_overlap(overlap)
	{
	}

	/** \brief return singularity type */
	singularity_type const &get_singularity_type(void) const
	{
		return m_sing_type;
	}

	/** \brief return overlapping state */
	element_overlapping const &get_overlap(void) const
	{
		return m_overlap;
	}

private:
	/** \brief the singularity type */
	singularity_type const m_sing_type;
	/** \brief the overlapping state */
	element_overlapping const m_overlap;
};


/** \brief function to determine the overlapping state of two elements
 * \tparam test_field_t the test field type
 * \tparam trial_field_t the trial field type
 * \param [in] test_field the test field
 * \param [in] trial_field the trial field
 * \return the element match object
 */
template <class test_field_t, class trial_field_t>
element_match element_match_eval(
	field_base<test_field_t> const &test_field,
	field_base<trial_field_t> const &trial_field)
{
	bool const face_match_possible = std::is_same<
		typename test_field_t::elem_t,
		typename trial_field_t::elem_t
	>::value;

	if (face_match_possible)
		if (test_field.get_elem().get_id() == trial_field.get_elem().get_id())
			return element_match(FACE_MATCH);

	/** \todo disable these cases for the collocational formalism */
/*
	element_overlapping overlap(test_field.get_elem().get_overlapping(trial_field.get_elem()));

	// check edge match
	if (overlap.get_num() > 1)
		return element_match(EDGE_MATCH, overlap);

	// check corner match
	if (overlap.get_num() == 1)
		return element_match(CORNER_MATCH, overlap);
*/
	return element_match(REGULAR);
}

#endif // ELEMENT_MATCH_HPP_INCLUDED
