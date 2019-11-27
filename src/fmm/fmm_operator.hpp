/**
 * @file fmm_operator.hpp 
 * @brief FMM operator types and tags 
 * @ingroup fmm_ops
 */

#ifndef FMM_OPERATOR_HPP_INCLUDED
#define FMM_OPERATOR_HPP_INCLUDED

#include <tuple>

namespace NiHu
{
namespace fmm
{

namespace op_tags
{
	struct m2m {};
	struct l2l {};
	struct m2l {};
	
	struct p2m {};
	struct p2l {};
	struct l2p {};
	struct m2p {};
	
	struct p2p {};
	
	namespace internal 
	{
		typedef std::tuple<m2m, l2l, m2l, p2m, p2l, l2p, m2p, p2p> op_tag_order;
		
		template <class T, class Tuple>
		struct index;

		template <class T, class... Types>
		struct index<T, std::tuple<T, Types...> > 
		{
			static const size_t value = 0;
		};

		template <class T, class U, class... Types>
		struct index<T, std::tuple<U, Types...> > 
		{
			static const size_t value = 1 + index<T, std::tuple<Types...>>::value;
		};
	}
	
	constexpr size_t num_tags()
	{
		return std::tuple_size<internal::op_tag_order>::value;
	}
	
	template <class Tag>
	constexpr size_t tag2idx(Tag const &tag = Tag())
	{
		return internal::index<Tag, internal::op_tag_order>::value;
	}
	
	template <size_t idx>
	struct idx2tag
	{
		typedef typename std::tuple_element<idx, internal::op_tag_order>::type type;
	};
}

/**
 * @todo These typedefs support previous syntax
 */
typedef op_tags::m2m m2m_tag;
typedef op_tags::l2l l2l_tag;
typedef op_tags::m2l m2l_tag;

typedef op_tags::m2p m2p_tag;
typedef op_tags::l2p l2p_tag;
typedef op_tags::p2l p2l_tag;
typedef op_tags::p2m p2m_tag;

typedef op_tags::p2p p2p_tag;

/** @brief Operator defining its tag type */
template <class FmmTag>
struct fmm_operator
{
	typedef FmmTag fmm_tag;
};

} // end of namespace fmm
} // end of namespace NiHu


#endif
