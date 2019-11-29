#ifndef NIHU_LOCAL_OPERATOR_HPP
#define NIHU_LOCAL_OPERATOR_HPP

#include <type_traits>

namespace NiHu
{
namespace fmm
{

template <class Operator>
struct is_local_operator
	: public std::false_type {};

}
}

#endif
