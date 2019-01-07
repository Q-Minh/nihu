/**
 * \file if.hpp 
 * \ingroup tmp
 * \brief Implementation of the \ref if_ metafunction
 */

#ifndef IF_HPP_INCLUDED
#define IF_HPP_INCLUDED

namespace tmp
{
	/**
	 * \brief IF control structure
	 * \tparam Choice a choice evaluated to a logical type
	 * \tparam T the type returned when Choice is true_type
	 * \tparam F the type returned when Choice is false_type
	 */
	template <class Cond, class T, class F>
	struct if_ : std::conditional<Cond::value, T, F> {};
}

#endif /* IF_HPP_INCLUDED */