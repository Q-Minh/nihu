/**
 * \file integral.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief Declaration of class Integral and its specialisations
 */
#ifndef INTEGRAL_HPP_INCLUDED
#define INTEGRAL_HPP_INCLUDED

#include "element.hpp"
#include <algorithm>

/**
 * \brief 
 */
template <class Domain, class Function>
class Integral
{
public:
	ret_t eval(Domain const &domain)
	{
		if (is_linear<Domain>::value)
			return std::accumulate(
				quad.begin(),
				quad.end(),
				[] (arg_t &q) { Kernel::eval(q.x) * q.w; }
			) * domain.get_jac();
		else
			return std::accumulate(
				quad.begin(),
				quad.end(),
				[] (arg_t &q) { Kernel::eval(q)  * q.w * domain.get_jac(q); }
			);
	}
};

#endif
