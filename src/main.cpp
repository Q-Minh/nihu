/**
 * \file main.cpp
 * \brief test file for numerical integration with quadrature
 */
 
#include "integrate.hpp"

#include <iostream>
#include <vector>

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
		apply (kernel_arg_type<myKernel>::type x)
	{
		return x*x;
	}
};


template <> struct product_type<double, double> { typedef double type; };

int main(void)
{
	std::vector<myQuad> quadrature;

	IntegFunctor<myKernel, myQuad>::res_t res = std::accumulate(quadrature.begin(), quadrature.end(), IntegFunctor<myKernel, myQuad>::res_t(), IntegFunctor<myKernel, myQuad>());

	std::cout << res << quadrature.size() << std::endl;

	return 0;
}
