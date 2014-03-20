#include "core/integral_operator.hpp"
#include "library/laplace_kernel.hpp"
#include <type_traits>
#include <iostream>

int main(void)
{
	// create an integral operator from a kernel
	auto intop = create_integral_operator(laplace_3d_SLP_kernel());
	auto k = intop.get_kernel();
	
	// create a couple integral operator from a set of kernels
	auto c_intop = create_integral_operator(
		laplace_3d_SLP_kernel(),
		laplace_3d_DLP_kernel());
	auto ck = c_intop.get_kernel();
	std::cout << std::boolalpha << is_couple<decltype(ck)>::value << std::endl;

	return 0;
}

