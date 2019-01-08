#include "../tmp/interval.hpp"

namespace tmp
{

/** \brief error terminating case of eval_interval */
template <>
int eval_interval<tmp::vector<> >(double)
{
	throw std::out_of_range("cannot determine interval value");
}

} // end of namespace tmp
