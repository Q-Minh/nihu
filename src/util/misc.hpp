/// \file misc.hpp
/// \brief miscellaneous functions and metafunctions
#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

namespace NiHu
{


/** \brief move selected elements from a range to an other
 * \tparam InputIterator the input iterator
 * \tparam OutputIterator the output iterator
 * \tparam Predicate the predicate
 * \param [in] first begin iterator of input range
 * \param [in] last end iterator of input range
 * \param [in] result begin iterator of output range
 * \param [in] pred the predicate
 */
template <class InputIterator, class OutputIterator, class Predicate>
InputIterator move_if(InputIterator first, InputIterator last, OutputIterator result, Predicate pred)
{
	InputIterator to = first;
	while (first != last)
	{
		if (pred(*first))
			*result++ = *first++;
		else
			*to++ = *first++;
	}
	return to;
}


} // namespace NiHu

#endif
