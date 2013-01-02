#include <iostream>

#include <vector>
#include <numeric>

template <class Quad>
struct quad_weight_type;

template <class Quad>
struct quad_location_type;

struct myQuad;
template <> struct quad_weight_type<myQuad> { typedef double type; };
template <> struct quad_location_type<myQuad> { typedef double type; };

struct myQuad {
	quad_location_type<myQuad>::type x;
	quad_weight_type<myQuad>::type w;
};



template <class Kernel>
struct kernel_result_type;

template <class Kernel>
struct kernel_arg_type;

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


template <class A, class B>
struct product_type;

template <> struct product_type<double, double> { typedef double type; };


template <class Kernel, class Quad>
struct IntegFunctor
{
	typedef typename product_type<
		typename kernel_result_type<Kernel>::type,
		typename quad_weight_type<Quad>::type
	>::type res_t;
	res_t operator () (res_t const &x, Quad const &q) 	{ return x + Kernel::apply(q.x)*q.w; }
};


int main(void)
{
	std::vector<myQuad> quadrature;

	IntegFunctor<myKernel, myQuad>::res_t res = std::accumulate(quadrature.begin(), quadrature.end(), IntegFunctor<myKernel, myQuad>::res_t(), IntegFunctor<myKernel, myQuad>());

	std::cout << res << quadrature.size() << std::endl;

	return 0;
}
