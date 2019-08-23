
#ifndef FMM_OPERATOR_COLLECTION_HPP_INCLUDED
#define FMM_OPERATOR_COLLECTION_HPP_INCLUDED

#include <type_traits>
#include <tuple>

namespace NiHu
{
namespace fmm
{

template <class ...Ops>
class fmm_operator_collection;

template <class ...Ops>
fmm_operator_collection<Ops...>
create_fmm_operator_collection(Ops &&... ops);

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

public:
	template <class FmmTag>
	struct op_type : op_type_impl<FmmTag, Ops...> {};

	fmm_operator_collection(Ops &&... ops)
		: m_ops(std::forward<Ops>(ops)...)
	{
	}

	template <typename FmmTag>
	auto get(FmmTag tag) const
	{
		return std::get<
			typename op_type<FmmTag>::type
		>(m_ops);
	}

private:
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

}
}

#endif
