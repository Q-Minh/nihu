/**
 * @file fmm_operator_collection.hpp 
 * @brief Implementation of @ref NiHu::fmm::fmm_operator_collection
 * @ingroup fmm_ops
 */

#ifndef FMM_OPERATOR_COLLECTION_HPP_INCLUDED
#define FMM_OPERATOR_COLLECTION_HPP_INCLUDED

#include <type_traits>
#include <tuple>

namespace NiHu
{
namespace fmm
{

// Forward declaration of class
template <class ...Ops>
class fmm_operator_collection;

// Forward declaration of creator function
template <class ...Ops>
fmm_operator_collection<Ops...>
create_fmm_operator_collection(Ops &&... ops);

/**
 * @brief Class representing a collection of FMM operators
 * @tparam ...Ops Contained operator types 
 * @details
 * The operator collection enables performing the same manipulations, such as
 * integration, indexing and precomputing on a group of FMM operators at the
 * same time.
 */ 
template <class ...Ops>
class fmm_operator_collection
{
private:
	template <class FmmTag, class First, class...Tail>
	struct op_type_impl : std::conditional<
		std::is_same<typename std::decay<First>::type::fmm_tag, FmmTag>::value,
		First,
		typename op_type_impl<FmmTag, Tail...>::type
	> {};

	template <class FmmTag, class Last>
	struct op_type_impl<FmmTag, Last> : std::conditional<
		std::is_same<typename std::decay<Last>::type::fmm_tag, FmmTag>::value,
		Last,
		void
	> {};

	/**
	 * @brief Initialize collection from tuple of operators 
	 * @details
	 * Only used internally for concatenation
	 */
	fmm_operator_collection(std::tuple<Ops...> && ops)
		: m_ops(std::forward<std::tuple<Ops...> >(ops))
	{
	}
	
public:
	/**
	 * @brief Metafunction for retreiving operator type for a given tag 
	 * @tparam FmmTag Operator type tag
	 * @details
	 * This metafunction is utilized for extracting the specified operator type
	 * from the collection.
	 */
	template <class FmmTag>
	struct op_type : op_type_impl<FmmTag, Ops...> {};

	/** 
	 * @brief Constructor 
	 * @param[in] ... Operators contained in the collection 
	 */
	fmm_operator_collection(Ops &&... ops)
		: m_ops(std::forward<Ops>(ops)...)
	{
	}

	/**
	 * @brief Retrieve operator from the collection
	 * @tparam FmmTag Operator type tag
	 * @param[in] tag Tag instance to identify operator type 
	 * @return Requested operator
	 * 
	 * @details 
	 * Static assertion guarantees that a compilation error occurs if the 
	 * requested operator is not found in the collection.
	 */
	template <typename FmmTag>
	auto get(FmmTag tag) const
	{
		typedef typename op_type<FmmTag>::type operator_t;
		
		// The required operator type must be present
		static_assert(!std::is_same<operator_t, void>::value,
					"The required operator type was not found in the collection");
		
		return std::get<operator_t>(m_ops);
	}
	
	/**
	 * @brief Concatenate (create union of) two operator collections
	 * @tparam ...OpsRhs Operators of the right hand side collection
	 * @return Concatenated collection
	 */
	template <class ...OpsRhs>
	auto operator|(fmm_operator_collection<OpsRhs...> &&rhs) const
	{
		return fmm_operator_collection(std::tuple_cat(m_ops, rhs.m_ops));
	}

	/**
	 * @brief Concatenate (create union of) a collection and an operator
	 * @tparam Operator Operator added to the collection
	 * @return Extended collection
	 */
	template <class Operator>
	auto operator|(Operator &&rhs) const
	{
		return this->operator|(create_fmm_operator_collection(std::forward<Operator>(rhs)));
	}
	
private:
	/**
	 * @brief Implementation of transform
	 * @tparam F Functor of transformation
	 * @tparam ...Is Index sequence for indexing the contained tuple
	 * @return Transformed collection
	 */
	template <class F, size_t... Is>
	auto transform_impl(F f, std::index_sequence<Is...>) {
		return create_fmm_operator_collection(
			f(std::get<Is>(m_ops))...
		);
	}

public:
	template <class F>
	auto transform(F f) {
		return transform_impl(
			f, std::make_index_sequence<sizeof...(Ops)>{});
	}

private:
	std::tuple<Ops...> m_ops;
};

template <class ...Ops>
fmm_operator_collection<Ops...>
create_fmm_operator_collection(Ops &&... ops)
{
	return fmm_operator_collection<Ops...>(std::forward<Ops>(ops)...);
}

} // end of namespace fmm
} // end of namespace NiHu

template <class Operator, class ...OpsColl>
auto operator|(Operator &&op, NiHu::fmm::fmm_operator_collection<OpsColl...> &&coll)
{
	return std::forward<NiHu::fmm::fmm_operator_collection<OpsColl...> >(coll) | std::forward<Operator>(op);
}


#endif /* FMM_OPERATOR_COLLECTION_HPP_INCLUDED */
