/**
* \file integral.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of class Integral and its specialisations
*/
#ifndef INTEGRAL_HPP_INCLUDED
#define INTEGRAL_HPP_INCLUDED

#include "elem_descriptor.hpp"
#include <numeric>

/**
* \brief 
*/
template <class ElemDescriptor, class Kernel>
class integral
{
public:
	typedef double ret_t;
	typedef Kernel kernel_t;

	template <class InputIterator>
	static ret_t eval(InputIterator begin, InputIterator end)
	{
		// initialise result to zero
		return std::accumulate(
			begin, end, ret_t(),
			[] (ret_t const &x, ElemDescriptor const &ed) { return x + kernel_t::eval(ed) * ed.get_jacobian(); }
		);
	}
};

#endif
