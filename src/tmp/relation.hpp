#ifndef RELATION_HPP_INCLUDED
#define RELATION_HPP_INCLUDED

namespace tmp
{
	/** \brief return true_type if first is less than second */
	template <class N, class M> struct less;

	/** \brief return true_type if first is greater than second */
	template <class N, class M> struct greater;


	/** \brief compute maximum of types */
	template <class Val, class...Args>
	struct max_ : max_<Val, typename max_<Args...>::type> {};

	/** \brief specialisation of max_ for the two parameter case */
	template <class Val1, class Val2>
	struct max_<Val1, Val2> : std::conditional<
		greater<Val1, Val2>::value,
		Val1, Val2
	> {};

	/** \brief compute minimum of types */
	template <class Val, class...Args>
	struct min_ : min_<Val, typename min_<Args...>::type> {};

	/** \brief specialisation of min_ for the two parameter case */
	template <class Val1, class Val2>
	struct min_<Val1, Val2> : std::conditional<
		less<Val1, Val2>::value,
		Val1, Val2
	> {};
}


#endif // RELATION_HPP_INCLUDED
