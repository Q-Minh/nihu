/**
 * \file kernel.hpp
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "descriptor.hpp"

#include <complex>
typedef std::complex<double> dcomplex;

/**
 * \brief base class of 3D Helmholtz kernels
 * \tparam NumElements number of scalars in the kernel
 */
template <unsigned NumElements>
class green_base
{
public:
	/** \brief template parameter as nested constant */
	static const unsigned num_elements = NumElements;

	/** \brief 3D location */
	typedef typename tria_1_elem::x_t x_t;
	/** \brief the scalar type of the kernel */
	typedef dcomplex scalar_t;
	/** \brief the result type of the kernel */
	typedef Eigen::Matrix<scalar_t, 1, num_elements> result_t;
	
	/** \brief set the wave number do a defined value */
	static void set_wave_number(dcomplex const &k)
	{
		green_base::k = k;
	}

	/** \brief set the source location to a defined value */
	static void set_x0(x_t const &x0)
	{
		green_base::x0 = x0;
	}
	
protected:
	/** \brief wave number */
	static dcomplex k;
	/** \brief source location */
	static x_t x0;
	/** \brief kernel result */
	static result_t result;
};

/** \brief definition of static member wave number */
template <unsigned NumElements>
dcomplex green_base<NumElements>::k;
/** \brief definition of static member source location */
template <unsigned NumElements>
typename green_base<NumElements>::x_t green_base<NumElements>::x0;
/** \brief definition of static member kernel result */
template <unsigned NumElements>
typename green_base<NumElements>::result_t green_base<NumElements>::result;


/**
 * \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$
 * \tparam NumElements number of scalars in the kernel
 */
class green_G_kernel : public green_base<1>
{
public:
	/** \brief kernel input type */
	typedef location<x_t> input_t;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		// source receiver distance
		double r = (input.get_x() - x0).norm();
		// complex kernel
		result[0] = std::exp(-dcomplex(0.0,1.0)*k*r) / r / 4.0 / M_PI;
		return result;
	}
};


/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$
 * \tparam NumElements number of scalars in the kernel
 */
class green_HG_kernel : public green_base<2>
{
public:
	/** \brief kernel input type */
	typedef location_with_normal<x_t> input_t;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		x_t rvec = input.get_x() - x0;

		double r = rvec.norm();
		result[0] = std::exp(-dcomplex(0.,1.)*k*r) / r / 4.0 / M_PI;
		// dr / dn
		double rdn = rvec.dot(input.get_normal()) / r;
		result[1] = result[0] * (-(1.0 + dcomplex(0.,1.)*k*r) / r) * rdn;
		return result;
	}
};

#endif

