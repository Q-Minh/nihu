#include "../library/normal_derivative_kernel.hpp"
#include "../library/laplace_base_kernel.hpp"
#include <iostream>

int main(void)
{
	typedef NiHu::space_2d<double> Space;
	typedef NiHu::laplace_helper<Space> Base;
	typedef NiHu::normal_derivative_kernel<Base, 0, 0> SLP;
	typedef NiHu::normal_derivative_kernel<Base, 0, 1> DLP;
	typedef NiHu::normal_derivative_kernel<Base, 1, 0> DLPt;
	typedef NiHu::normal_derivative_kernel<Base, 1, 1> HSP;
	
	Base b;
	SLP slp(b);
	DLP dlp(b);
	DLPt dlpt(b);
	HSP hsp(b);
	
	SLP::x_t x, y, nx, ny;
	x << 0.0, 0.0;
	y << 1.0, 0.0;
	nx << 0.0, 1.0;
	ny << 1.0, 0.0;
	
	std::cout << "SLP: " << slp(x, y) << std::endl;
	std::cout << "DLP: " << dlp(x, y, ny) << std::endl;
	std::cout << "DLPt: " << dlpt(x, y, nx) << std::endl;
	std::cout << "HSP: " << hsp(x, y, nx, ny) << std::endl;
	
	return 0;
}
