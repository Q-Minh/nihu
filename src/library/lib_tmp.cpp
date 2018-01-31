#include "../tmp/interval.hpp"

namespace NiHu
{

/** \brief error terminating case of eval_interval */
template <>
int eval_interval<tmp::vector<> >(double)
{
	throw std::out_of_range("cannot determine interval value");
}

}
