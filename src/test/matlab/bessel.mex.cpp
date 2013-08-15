#include "util/eigen_utils.hpp"
#include "util/math_functions.hpp"
#include <iostream>

int main(void)
{
	for (int i = 0; i < 10; ++i)
	{
		auto z = std::complex<double>(0.0,1.0)*(1.0*i);
		auto H0 = bessel::H<0>(z);
		auto H1 = bessel::H<1>(z);

		std::cout << z << '\t' << H0 << '\t' << H1 << std::endl;
	}
}
