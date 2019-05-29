/** \file type2tag.hpp 
 * \brief assign a tag to a type
 */
#ifndef FMM_TYPE_TO_TAG_HPP_INCLUDED
#define FMM_TYPE_TO_TAG_HPP_INCLUDED

namespace NiHu
{
namespace fmm
{

/** \brief metafunction to assign a tag to a type
 * \tparam Type the input type
 */
template<class Type>
struct type2tag {};

/** \brief metafunction to recover the type from a tag
 * \tparam Tag the input tag
 */
template <class Tag>
struct tag2type;

template <class Type>
struct tag2type<type2tag<Type> >
{
	typedef Type type;
};

} // end of namespace fmm
} // namespace NiHu

#endif // FMM_TYPE_TO_TAG_HPP_INCLUDED
