#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "descriptor.hpp"

#include <complex>
typedef std::complex<double> dcomplex;

class green_G_kernel
{
public:
	typedef location<typename tria_1_elem::x_t> input_t;
	typedef dcomplex scalar_t;
	static const unsigned num_elements = 1;
	typedef Eigen::Matrix<scalar_t, 1, num_elements> result_t;

	typedef input_t::x_t x_t;

	static void set_wave_number(dcomplex const &k)
	{
		green_G_kernel::k = k;
	}

	static void set_x0(x_t const &x0)
	{
		green_G_kernel::x0 = x0;
	}

	static result_t const &eval (input_t const &input)
	{
		double r = (input.get_x() - x0).norm();
		result[0] = std::exp(-dcomplex(0.0,1.0)*k*r) / r / 4.0 / M_PI;
		return result;
	}

protected:
	static dcomplex k;
	static x_t x0;
	static result_t result;
};

dcomplex green_G_kernel::k;
green_G_kernel::x_t green_G_kernel::x0;
green_G_kernel::result_t green_G_kernel::result;


class green_HG_kernel
{
public:
	typedef location_with_normal<typename tria_1_elem::x_t> input_t;
	typedef dcomplex scalar_t;
	static const unsigned num_elements = 2;
	typedef Eigen::Matrix<scalar_t, 1, num_elements> result_t;

	typedef input_t::x_t x_t;

	static void set_wave_number(dcomplex const &k)
	{
		green_HG_kernel::k = k;
	}

	static void set_x0(x_t const &x0)
	{
		green_HG_kernel::x0 = x0;
	}

	static result_t const &eval (input_t const &input)
	{
		x_t rvec = input.get_x() - x0;
		x_t normal = input.get_normal();
		normal /= normal.norm();

		double r = rvec.norm();
		result[0] = std::exp(-dcomplex(0.,1.)*k*r) / r / 4.0 / M_PI;
		result[1] = result[0] * (-(1.0 + dcomplex(0.,1.)*k*r) / r) * rvec.dot(normal) / r;
		return result;
	}

protected:
	static dcomplex k;
	static x_t x0;
	static result_t result;
};

dcomplex green_HG_kernel::k;
green_HG_kernel::x_t green_HG_kernel::x0;
green_HG_kernel::result_t green_HG_kernel::result;

#endif

