/**
 * \file kernel.hpp
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "descriptor.hpp"

#include <complex>
/** \brief double complex type */
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
		green_base::m_k = k;
	}

	/** \brief set the source location to a defined value */
	static void set_x0(x_t const &x0)
	{
		green_base::m_x0 = x0;
	}

	/**
	 * \brief associate relative distance with a required polynomial degree
	 */
	struct kernel_precision {
		double rel_distance;	/**< \brief the relative distance from the kernel's source point */
		unsigned degree;		/**< \brief polynomial degree required for accurate integration */
	};

protected:
	/** \brief wave number */
	static dcomplex m_k;
	/** \brief source location */
	static x_t m_x0;
	/** \brief kernel result */
	static result_t m_result;
};

/** \brief definition of static member wave number */
template <unsigned NumElements>
dcomplex green_base<NumElements>::m_k;
/** \brief definition of static member source location */
template <unsigned NumElements>
typename green_base<NumElements>::x_t green_base<NumElements>::m_x0;
/** \brief definition of static member kernel result */
template <unsigned NumElements>
typename green_base<NumElements>::result_t green_base<NumElements>::m_result;


/**
 * \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$
 * \tparam NumElements number of scalars in the kernel
 */
class green_G_kernel : public green_base<1>
{
public:
	typedef green_base<1> base;	/**< \brief the base class' type */

	/** \brief kernel input type */
	typedef location<x_t> input_t;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		// source receiver distance
		double r = (input.get_x() - m_x0).norm();
		// complex kernel
		m_result[0] = std::exp(-dcomplex(0.0,1.0)*m_k*r) / r / (4.0 * M_PI);
		return m_result;
	}

	/** \brief estimate kernel complexity for a given input
	 * \param input
	 * \return plynomial order of the kernel
	 */
	static unsigned estimate_complexity(input_t const &input)
	{
		double rel_distance = ((input.get_x() - m_x0).norm()) / sqrt(input.get_jacobian());
		for (unsigned i = 0; i < 6; ++i)
			if (limits[i].rel_distance < rel_distance)
				return limits[i].degree;
		return 5;
	}

protected:
	static const base::kernel_precision limits[];	/**< \brief array of distance limits */
};

const green_G_kernel::base::kernel_precision green_G_kernel::limits[6] = {
	{9.2, 1},
	{1.6, 3},
	{1.0, 5},
	{0.8, 7},
	{0.7, 9},
	{0.6, 11}
	};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$
 * \tparam NumElements number of scalars in the kernel
 */
class green_HG_kernel : public green_base<2>
{
public:
	typedef green_base<2> base;	/**< \brief base class' type */

	/** \brief kernel input type */
	typedef location_with_normal<x_t> input_t;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		x_t rvec = input.get_x() - m_x0;
		double r2 = rvec.squaredNorm();
		double r = sqrt(r2);
		dcomplex ikr(dcomplex(0.,1.)*m_k*r);

		m_result[0] = std::exp(-ikr) / r / (4.0 * M_PI);
		double rdn = rvec.dot(input.get_normal());
		m_result[1] = m_result[0] * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

	/** \brief estimate kernel complexity for a given input
	 * \param input
	 * \return plynomial order of the kernel
	 */
	static unsigned estimate_complexity(input_t const &input)
	{
		double rel_distance = ((input.get_x() - m_x0).norm()) / sqrt(input.get_jacobian());
		for (unsigned i = 0; i < 3; ++i)
			if (limits[i].rel_distance < rel_distance)
				return limits[i].degree;
		return 7;
	}

protected:
	static const base::kernel_precision limits[];	/**< \brief array of distance limits */
};

const green_HG_kernel::base::kernel_precision green_HG_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{1.5, 7}
	};

#endif

