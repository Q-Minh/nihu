#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "descriptor.hpp"

#include <complex>
typedef std::complex<double> dcomplex;

class green_kernel
{
public:
	typedef location<Eigen::Matrix<double, 1, 3> > input_t;
	typedef dcomplex result_t;

	typedef input_t::x_t x_t;

	static void set_wave_number(dcomplex const &k)
	{
		green_kernel::k = k;
	}

	static void set_x0(x_t const &x0)
	{
		green_kernel::x0 = x0;
	}

	static result_t eval (input_t const &input)
	{
		double r = (input.get_x() - x0).norm();
		return std::exp(-dcomplex(0.0,1.0)*k*r) / r;
	}

protected:
	static dcomplex k;
	static x_t x0;
};

dcomplex green_kernel::k;
Eigen::Matrix<double, 1, 3> green_kernel::x0;

#endif

