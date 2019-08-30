/** 
 * \file type2tag.hpp 
 * \brief Assign a tag to a type
 * \ingroup util
 */
#ifndef TYPE2TAG_HPP_INCLUDED
#define TYPE2TAG_HPP_INCLUDED

namespace NiHu
{

/** 
 * \brief Netafunction assigning a tag to a type
 * \tparam Type Input type
 */
template<class Type>
struct type2tag 
{
	// self-returning metafunction
	typedef type2tag type;
};

/** 
 * \brief Metafunction recovering the type from a tag
 * \tparam Tag Input tag
 */
template <class Tag>
struct tag2type;

template <class Type>
struct tag2type<type2tag<Type> >
{
	typedef Type type;
};

} // end of namespace NiHu

#endif /* TYPE2TAG_HPP_INCLUDED */
