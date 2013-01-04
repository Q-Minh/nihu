/**
 * \file main.cpp
 * \brief test file for numerical integration with quadrature
 */
 
#include "integrate.hpp"

#include <iostream>

struct myQuad;
template <> struct quad_weight_type<myQuad> { typedef double type; };
template <> struct quad_location_type<myQuad> { typedef double type; };

struct myQuad {
	quad_location_type<myQuad>::type x;
	quad_weight_type<myQuad>::type w;
};


struct myKernel;
template<> struct kernel_result_type<myKernel> { typedef double type; };
template<> struct kernel_arg_type<myKernel> { typedef double type; };

struct myKernel
{
	static kernel_result_type<myKernel>::type
		eval (kernel_arg_type<myKernel>::type x)
	{
		return x*x;
	}
};


int main(void)
{
	myQuad q[] = {
		{-1.0, 1.0/3.0},
		{0.0, 1.0/3.0},
		{1.0, 1.0/3.0}
	};

	std::cout << Integral<myQuad*, myKernel>::eval(q, q+3) << std::endl;

	return 0;
}

