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

// forward declaration
template <class Derived>
struct kernel_traits;

/**
 * \brief CRTP base class of all BEM kernels
 * \tparam the CRTP derived class
 */
template <class Derived>
class kernel_base
{
public:
	/** \brief type of the location */
	typedef typename kernel_traits<Derived>::x_t x_t;
	/** \brief type of the kernel input */
	typedef typename kernel_traits<Derived>::input_t input_t;
	/** \brief type of the kernel's scalar */
	typedef typename kernel_traits<Derived>::scalar_t scalar_t;
	/** \brief type of the kernel's result */
	typedef typename kernel_traits<Derived>::result_t result_t;

	static void reset(void)
	{
		m_num_evaluations = 0;
	}

	static unsigned get_num_evaluations(void)
	{
		return m_num_evaluations;
	}

	/**
	 * \brief set the source location number to a defined value
	 * \param [in] x0 source location
	 */
	static void set_x0(x_t const &x0)
	{
		m_x0 = x0;
	}

	/**
	 * \brief evaluate kernel at a given receiver position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x
	 */
	static result_t const &eval (input_t const &x)
	{
		m_num_evaluations++;
		return Derived::eval(x);
	}

	/**
	 * \brief evaluate kernel at a given source and receiver position
	 * \param [in] x0 source position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x0 and x
	 */
	static result_t const &eval (input_t const &x0, input_t const &x)
	{
		set_x0(x0.get_x());
		return eval(x);
	}

protected:
	static x_t m_x0;			/**< \brief source location */
	static result_t m_result;	/**< \brief kernel result */
	static unsigned m_num_evaluations;	/**< \brief number of kernel evaluations */
};

/**< \brief static member source location */
template <class Derived>
typename kernel_base<Derived>::x_t kernel_base<Derived>::m_x0;
/**< \brief static member kernel result */
template <class Derived>
typename kernel_base<Derived>::result_t kernel_base<Derived>::m_result;
/**< \brief static member kernel result */
template <class Derived>
unsigned kernel_base<Derived>::m_num_evaluations = 0;


/**
 * \brief CRTP base class of 3D Helmholtz kernels and its derivatives
 * \tparam Derived type of the derived class in the CRTP scheme
 */
template <class Derived>
class helmholtz_base : public kernel_base<Derived>
{
public:
	typedef kernel_base<Derived> base_t;
	using x_t = typename base_t::x_t;
	using input_t = typename base_t::input_t;
	using scalar_t = typename base_t::scalar_t;
	using result_t = typename base_t::result_t;

	using base_t::set_x0;
	using base_t::eval;
	using base_t::m_x0;

	/** \brief associate relative distance with a required polynomial degree */
	struct kernel_precision
	{
		double rel_distance;	/**< \brief relative distance from the source point */
		unsigned degree;		/**< \brief polynomial degree required for integration */
	};

	/**
	 * \brief set the wave number to a defined value
	 * \param [in] k wave number
	 */
	static void set_wave_number(dcomplex const &k)
	{
		helmholtz_base::m_k = k;
	}

	/**
	 * \brief determine kernel's polynomial complexity at a given receiver position
	 * \param [in] x receiver position
	 * \return polynomial degree needed for accurate integration
	 */
	static unsigned estimate_complexity(input_t const &x)
	{
		double rel_distance = ((x.get_x() - m_x0).norm()) / sqrt(x.get_jacobian());
		unsigned i;
		for (i = 0; rel_distance < Derived::limits[i].rel_distance; ++i);
		return Derived::limits[i].degree;
	}

	/**
	 * \brief determine kernel's polynomial complexity at a given receiver position
	 * \param [in] x0 source position
	 * \param [in] x receiver position
	 * \return polynomial degree needed for accurate integration
	 */
	static unsigned estimate_complexity(input_t const &x0, input_t const &x)
	{
		set_x0(x0.get_x());
		return estimate_complexity(x);
	}

protected:
	static dcomplex m_k;		/**< \brief wave number */
};

/**< \brief static member wave number */
template <class Derived>
dcomplex helmholtz_base<Derived>::m_k;


// forward declaration
class helmholtz_G_kernel;

/** \brief traits of the helmholtz G kernel */
template<>
struct kernel_traits<helmholtz_G_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location<x_t> input_t;
	typedef dcomplex result_t;
	typedef gauss_family_tag quadrature_family_t;
};

/**
 * \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$
 */
class helmholtz_G_kernel : public helmholtz_base<helmholtz_G_kernel>
{
	friend class helmholtz_base<helmholtz_G_kernel>;
public:
	typedef helmholtz_base<helmholtz_G_kernel> base_t;	/**< \brief the base class' type */

	using base_t::eval;

	/**
	 * \brief evaluate kernel at a given receiver position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x
	 */
	static result_t const &eval(input_t const &x)
	{
		// source receiver distance
		double r = (x.get_x() - m_x0).norm();
		// complex kernel
		m_result = std::exp(-dcomplex(0.0,1.0)*m_k*r) / r / (4.0 * M_PI);
		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const helmholtz_G_kernel::kernel_precision helmholtz_G_kernel::limits[] = {
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

class helmholtz_H_kernel;

/** \brief traits of the Helmholtz H kernel */
template<>
struct kernel_traits<helmholtz_H_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location_with_normal<x_t> input_t;
	typedef dcomplex result_t;
	typedef gauss_family_tag quadrature_family_t;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 */
class helmholtz_H_kernel : public helmholtz_base<helmholtz_H_kernel>
{
	friend class helmholtz_base<helmholtz_H_kernel>;
public:
	typedef helmholtz_base<helmholtz_H_kernel> base_t;	/**< \brief the base class' type */

	using base_t::eval; // needed so that the class sees the two argument eval function too

	/**
	 * \brief evaluate kernel at a given receiver position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x
	 */
	static result_t const &eval(input_t const &x)
	{
		x_t rvec = x.get_x() - m_x0;
		double r2 = rvec.squaredNorm();
		double r = sqrt(r2);
		dcomplex ikr(dcomplex(0.,1.)*m_k*r);

		m_result = std::exp(-ikr) / r / (4.0 * M_PI);
		double rdn = rvec.dot(x.get_normal());
		m_result = m_result * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const helmholtz_H_kernel::kernel_precision helmholtz_H_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};


class helmholtz_HG_kernel;

/** \brief traits of the double helmholtz kernel */
template<>
struct kernel_traits<helmholtz_HG_kernel>
{
	typedef tria_1_elem::x_t x_t;
	typedef dcomplex scalar_t;
	typedef location_with_normal<x_t> input_t;
	typedef couple<dcomplex> result_t;
	typedef gauss_family_tag quadrature_family_t;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$
 */
class helmholtz_HG_kernel : public helmholtz_base<helmholtz_HG_kernel>
{
	friend class helmholtz_base<helmholtz_HG_kernel>;
public:
	typedef helmholtz_base<helmholtz_HG_kernel> base_t;	/**< \brief base class' type */

//	using base_t::eval;

	/**
	 * \brief evaluate kernel at a given receiver position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x
	 */
	static result_t const &eval(input_t const &x)
	{
		x_t rvec = x.get_x() - m_x0;
		double r2 = rvec.squaredNorm();
		double r = sqrt(r2);
		dcomplex ikr(dcomplex(0.,1.)*m_k*r);

		m_result.first() = std::exp(-ikr) / r / (4.0 * M_PI);
		double rdn = rvec.dot(x.get_normal());
		m_result.second() = m_result.first() * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

protected:
	static const kernel_precision limits[];	/**< \brief array of distance limits */
};

const helmholtz_HG_kernel::kernel_precision helmholtz_HG_kernel::limits[]  = {
	{5.0, 2},
	{2.0, 5},
	{0.0, 7}
	};

#endif // KERNEL_HPP_INCLUDED

