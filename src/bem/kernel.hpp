/**
 * \file kernel.hpp
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "kernel_input.hpp"
#include "couple.hpp"

#include <complex>
/** \brief double complex type */
typedef std::complex<double> dcomplex;

template <class Derived>
struct green_kernel_traits;

/**
 * \brief base class of 3D Helmholtz kernels
 */
template <class Derived>
class green_base
{
public:
	typedef typename green_kernel_traits<Derived>::x_t x_t;
	typedef typename green_kernel_traits<Derived>::input_t input_t;
	typedef typename green_kernel_traits<Derived>::scalar_t scalar_t;
	typedef typename green_kernel_traits<Derived>::result_t result_t;

	/** \brief set the wave number to a defined value */
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

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input1)
	{
		return Derived::eval(input1);
	}

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input1, input_t const &input2)
	{
		set_x0(input1.get_x());
		return eval(input2);
	}

	static unsigned estimate_complexity(input_t const &input)
	{
		double rel_distance = ((input.get_x() - m_x0).norm()) / sqrt(input.get_jacobian());
		unsigned i;
		for (i = 0; rel_distance < Derived::limits[i].rel_distance; ++i);
		return Derived::limits[i].degree;
	}

	static unsigned estimate_complexity(input_t const &input1, input_t const &input2)
	{
		set_x0(input1.get_x());
		return estimate_complexity(input2);
	}

protected:
	/** \brief wave number */
	static dcomplex m_k;
	/** \brief source location */
	static x_t m_x0;
	/** \brief kernel result */
	static result_t m_result;
};

/** \brief definition of static member wave number */
template <class Derived>
dcomplex green_base<Derived>::m_k;
/** \brief definition of static member source location */
template <class Derived>
typename green_base<Derived>::x_t green_base<Derived>::m_x0;
/** \brief definition of static member kernel result */
template <class Derived>
typename green_base<Derived>::result_t green_base<Derived>::m_result;


class green_G_kernel;

template<>
struct green_kernel_traits<green_G_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location<x_t> input_t;
	typedef dcomplex result_t;
};

/**
 * \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$
 */
class green_G_kernel : public green_base<green_G_kernel>
{
	friend class green_base<green_G_kernel>;
public:
	typedef green_base<green_G_kernel> base_t;	/**< \brief the base class' type */

	using base_t::eval;
	using base_t::estimate_complexity;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		// source receiver distance
		double r = (input.get_x() - m_x0).norm();
		// complex kernel
		m_result = std::exp(-dcomplex(0.0,1.0)*m_k*r) / r / (4.0 * M_PI);
		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const green_G_kernel::kernel_precision green_G_kernel::limits[] = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
/*
	{9.2, 1},
	{1.6, 3},
	{1.0, 5},
	{0.8, 7},
	{0.7, 9},
	{0.6, 11}
*/
	};

class green_H_kernel;

template<>
struct green_kernel_traits<green_H_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location_with_normal<x_t> input_t;
	typedef dcomplex result_t;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 */
class green_H_kernel : public green_base<green_H_kernel>
{
	friend class green_base<green_H_kernel>;
public:
	typedef green_base<green_H_kernel> base_t;	/**< \brief the base class' type */

	using base_t::eval;
	using base_t::estimate_complexity;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		x_t rvec = input.get_x() - m_x0;
		double r2 = rvec.squaredNorm();
		double r = sqrt(r2);
		dcomplex ikr(dcomplex(0.,1.)*m_k*r);

		m_result = std::exp(-ikr) / r / (4.0 * M_PI);
		double rdn = rvec.dot(input.get_normal());
		m_result = m_result * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const green_H_kernel::kernel_precision green_H_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};


class green_HG_kernel;

template<>
struct green_kernel_traits<green_HG_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location_with_normal<x_t> input_t;
	typedef couple<dcomplex> result_t;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$
 */
class green_HG_kernel : public green_base<green_HG_kernel>
{
	friend class green_base<green_HG_kernel>;
public:
	typedef green_base<green_HG_kernel> base_t;	/**< \brief base class' type */

	using base_t::eval;
	using base_t::estimate_complexity;

	/** \brief evaluate the kernel for a given input */
	static result_t const &eval (input_t const &input)
	{
		x_t rvec = input.get_x() - m_x0;
		double r2 = rvec.squaredNorm();
		double r = sqrt(r2);
		dcomplex ikr(dcomplex(0.,1.)*m_k*r);

		m_result.first() = std::exp(-ikr) / r / (4.0 * M_PI);
		double rdn = rvec.dot(input.get_normal());
		m_result.second() = m_result.first() * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const green_HG_kernel::kernel_precision green_HG_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};

#endif

