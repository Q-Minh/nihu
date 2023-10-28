#include "nihu/core/integral_operator.hpp"
#include "nihu/library/laplace_kernel.hpp"
#include <type_traits>
#include <iostream>

int main(void)
{
	// create an integral operator from a kernel
	auto intop = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
	auto k = intop.get_kernel();
	
	(void)k;
	
	return 0;
}
