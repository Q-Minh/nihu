/**
 * \file kernel.hpp
 * \brief implementation of various kernels
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
#ifndef KERNEL_HPP_INCLUDED
#define KERNEL_HPP_INCLUDED

#include "kernel_input.hpp"
#include "couple.hpp"

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
	typedef kernel_traits<Derived> traits_t;
	
	/** \brief type of the kernel input */
	typedef typename traits_t::test_input_t test_input_t;
	typedef typename traits_t::trial_input_t trial_input_t;

	static_assert(std::is_same<typename test_input_t::space_t, typename trial_input_t::space_t>::value,
		"The test and trial kernel inputs must define the same coordinate space");

	/** \brief type of the space */
	typedef typename test_input_t::space_t space_t;

	/** \brief type of the location */
	typedef typename space_t::location_t x_t;

	/** \brief type of the scalar */
	typedef typename space_t::scalar_t scalar_t;

	/** \brief type of the kernel's result */
	typedef typename traits_t::result_t result_t;
	/** \brief quadrature family tag */
	typedef typename traits_t::quadrature_family_t quadrature_family_t;
	/** \brief shows if kernel is symmetric */
	static bool const is_symmetric = traits_t::is_symmetric;
	/** \brief the singularity order (r^-order) */
	static unsigned const singularity_order = traits_t::singularity_order;
	/** \brief the singular quadrature */
	static unsigned const singular_quadrature_order = traits_t::singular_quadrature_order;

private:
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	class kernel_bind
	{
	public:
		kernel_bind(Derived &kernel, test_input_t const &x)
			: m_kernel(kernel), m_test_input(x)
		{
		}

		result_t eval(trial_input_t const &y) const
		{
			return m_kernel.eval(m_test_input, y);
		}

		unsigned estimate_complexity(trial_input_t const &y, scalar_t const &reference_size) const
		{
			return m_kernel.estimate_complexity(m_test_input, y, reference_size);
		}

	private:
		Derived &m_kernel;
		test_input_t const &m_test_input; 
	};

	kernel_base() :
		m_num_evaluations(0)
	{
	}

	kernel_bind bind(test_input_t const &x)
	{
		return kernel_bind(derived(), x);
	}

	/**
	 * \brief return number of kernel evaluations
	 * \return number of kernel evaluations
	 */
	long long unsigned get_num_evaluations(void) const
	{
		return m_num_evaluations;
	}

	/**
	 * \brief evaluate kernel at a given source and receiver position
	 * \param [in] x0 source position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x0 and x
	 */
	result_t eval(test_input_t const &x, trial_input_t const &y)
	{
		m_num_evaluations++;
		return derived()(x, y);
	}

	/**
	 * \brief determine kernel's polynomial complexity at a given receiver position
	 * \param [in] x receiver position
	 * \return polynomial degree needed for accurate integration
	 */
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t const &reference_size) const
	{
		return derived().estimate_complexity(x, y, reference_size);
	}

protected:
	long long unsigned m_num_evaluations;	/**< \brief number of kernel evaluations */
	result_t m_result;	/**< \brief kernel result */
};


#include "gaussian_quadrature.hpp"


// forward declaration
class unit_kernel;

/** \brief traits of the 3D unit kernel */
template<>
struct kernel_traits<unit_kernel>
{
	/** \brief test input type */
	typedef location<space_3d> test_input_t;
	/** \brief trial input type */
	typedef location<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef double result_t;
	/** \brief quadrature family tag */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief shows if kernel is symmetric */
	static bool const is_symmetric = true;
	/** \brief the singularity order (r^-order) */
	static unsigned const singularity_order = 0;
	/** \brief the singular quadrature order */
	static unsigned const singular_quadrature_order = 0;
};

/**
 * \brief 3D unit kernel
 */
class unit_kernel : public kernel_base<unit_kernel>
{
public:
	typedef kernel_base<unit_kernel> base_t;
	
	/**
	 * \brief evaluate kernel at a given receiver position
	 * \param [in] x receiver position
	 * \return kernel evaluated in x
	 */
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return 1.0;
	}
	
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t const &s) const
	{
		return 0;
	}
};



// forward declaration
class helmholtz_G_kernel;

/** \brief traits of the helmholtz G kernel */
template<>
struct kernel_traits<helmholtz_G_kernel>
{
	typedef location<space_3d> test_input_t;
	typedef location<space_3d> trial_input_t;
	typedef std::complex<space_3d::scalar_t> result_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	static unsigned const singularity_order = 1;
	static unsigned const singular_quadrature_order = 7;
};

/**
 * \brief 3D Helmholtz kernel \f$\exp(-ikr)/4\pi r\f$
 */
class helmholtz_G_kernel : public kernel_base<helmholtz_G_kernel>
{
public:
	void set_wave_number(std::complex<scalar_t> const &k)
	{
		m_k = k;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		scalar_t r = (y.get_x() - x.get_x()).norm();
		return std::exp(-std::complex<scalar_t>(0.0,1.0)*m_k*r) / r / (4.0 * M_PI);
	}
	
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t size) const
	{
		/** \todo hard coding of this 3 is sick */
		return 3;
	}

private:
	std::complex<scalar_t> m_k;
};



class helmholtz_H_kernel;

/** \brief traits of the Helmholtz H kernel */
template<>
struct kernel_traits<helmholtz_H_kernel>
{
	typedef location<space_3d> test_input_t;
	typedef location_with_normal<space_3d> trial_input_t;
	typedef std::complex<space_3d::scalar_t> result_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	static unsigned const singularity_order = 1;
	static unsigned const singular_quadrature_order = 7;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left(-(1+ikr)/r\right) \cdot dr/dn \f$
 */
class helmholtz_H_kernel : public kernel_base<helmholtz_H_kernel>
{
public:
	typedef kernel_base<helmholtz_H_kernel> base_t;	/**< \brief the base class' type */

	void set_wave_number(std::complex<scalar_t> const &k)
	{
		m_k = k;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		x_t rvec = y.get_x() - x.get_x();
		scalar_t r2 = rvec.squaredNorm();
		scalar_t r = sqrt(r2);
		std::complex<scalar_t> ikr(std::complex<scalar_t>(0.,1.)*m_k*r);

		result_t m_result = std::exp(-ikr) / r / (4.0 * M_PI);
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result *= (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}

	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t size) const
	{
		/** \todo hard coding of this 5 is sick */
		return 5;
	}

private:
	std::complex<scalar_t> m_k;
};


class helmholtz_GH_kernel;

/** \brief traits of the double helmholtz kernel */
template<>
struct kernel_traits<helmholtz_GH_kernel>
{
	typedef location<space_3d> test_input_t;
	typedef location_with_normal<space_3d> trial_input_t;
	typedef couple<std::complex<space_3d::scalar_t> > result_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	static unsigned const singularity_order = 1;
	static unsigned const singular_quadrature_order = 7;
};

/**
 * \brief 3D Helmholtz kernel \f$ \exp(-ikr)/4\pi r \left\{1, -(1+ikr)/r \cdot dr/dn\right\} \f$
 */
class helmholtz_GH_kernel : public kernel_base<helmholtz_GH_kernel>
{
public:
	typedef kernel_base<helmholtz_GH_kernel> base_t;	/**< \brief base class' type */

	void set_wave_number(std::complex<scalar_t> const &k)
	{
		m_k = k;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		x_t rvec = y.get_x() - x.get_x();
		scalar_t r2 = rvec.squaredNorm();
		scalar_t r = sqrt(r2);
		std::complex<scalar_t> ikr(std::complex<scalar_t>(0.,1.)*m_k*r);

		result_t m_result(std::exp(-ikr) / r / (4.0 * M_PI));
		scalar_t rdn = rvec.dot(y.get_unit_normal());
		m_result.second() = m_result.first() * (1.0 + ikr) * (-rdn / r2);

		return m_result;
	}
	
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t size) const
	{
		/** \todo hard coding of this 5 is sick */
		return 5;
	}

private:
	std::complex<scalar_t> m_k;
};

#endif // KERNEL_HPP_INCLUDED

