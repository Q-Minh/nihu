/**
 * \file poisson_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Poisson equation \f$\nabla^2 p = 0\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef POISSON_KERNEL_HPP_INCLUDED
#define POISSON_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "location_normal.hpp"
#include "reciprocal_distance_kernel.hpp"


template <class space>
struct distance_vector_brick
{
	template <class wall>
	class brick : public wall
	{
	public:
		template <class test_input_t, class trial_input_t>
		brick(test_input_t const &test_input, trial_input_t const &trial_input) :
			wall(test_input, trial_input),
			m_distance_vector(trial_input.get_x()-test_input.get_x())
		{
		}
		
		typename space::location_t const &
			get_distance_vector(void) const
		{
			return m_distance_vector;
		}
		
	private:
		typename space::location_t m_distance_vector;
	};
};


template <class scalar>
struct distance_brick
{
	template <class wall>
	class brick : public wall
	{
	public:
		template <class test_input_t, class trial_input_t>
		brick(test_input_t const &test_input, trial_input_t const &trial_input) :
			wall(test_input, trial_input),
			m_distance(wall::get_distance_vector().norm())
		{
		}
		
		scalar const & get_distance(void) const
		{
			return m_distance;
		}
		
	private:
		scalar m_distance;
	};
};


template <class scalar>
struct poisson_g_brick
{
	template <class wall>
	class brick : public wall
	{
	public:
		template <class test_input_t, class trial_input_t>
		brick(test_input_t const &test_input, trial_input_t const &trial_input) :
			wall(test_input, trial_input),
			m_poisson_g(1.0 / wall::get_distance() / (4.0 * M_PI))
		{
		}
		
		scalar const & get_poisson_g(void) const
		{
			return m_poisson_g;
		}
		
		scalar const & get_result(void) const
		{
			return m_poisson_g;
		}
		
	private:
		scalar m_poisson_g;
	};
};


template <class space>
struct poisson_g_wall : build<
	distance_vector_brick<space>,
	distance_brick<typename space::scalar_t>,
	poisson_g_brick<typename space::scalar_t>
> {};


// forward declaration
class poisson_G_kernel;

/** \brief traits of the poisson G kernel */
template<>
struct kernel_traits<poisson_G_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef typename poisson_g_wall<space_3d>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$1/4\pi r\f$ */
class poisson_G_kernel :
	public kernel_base<poisson_G_kernel>,
	public reciprocal_distance_kernel<poisson_G_kernel>
{
public:
	/** \brief the crtp base's type */
	typedef kernel_base<poisson_G_kernel> base_t;
	/** \brief type of the test input */
	typedef base_t::test_input_t test_input_t;
	/** \brief type of the trial input */
	typedef base_t::trial_input_t trial_input_t;
	/** \brief type of the scalar of the inputs */
	typedef base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	template <class test_input_t, class trial_input_t>
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		// instantiate output
		typedef typename kernel_traits<poisson_G_kernel>::output_t output_t;
		output_t output(x, y);
		// select result from output
		return output.get_result();
	}
};


template <class scalar>
struct poisson_g_brick_immediate
{
	template <class wall>
	class brick : public wall
	{
	public:
		template <class test_input_t, class trial_input_t>
		brick(test_input_t const &test_input, trial_input_t const &trial_input) :
			wall(test_input, trial_input)
		{
		}
		
		scalar const & get_poisson_g(void) const
		{
			return m_poisson_g;
		}
		
		scalar const & get_result(void) const
		{
			return m_poisson_g;
		}
		
	private:
		scalar m_poisson_g;
	};
};


template <class space>
struct poisson_g_wall_immediate : build<
	poisson_g_brick_immediate<typename space::scalar_t>
> {};


class poisson_G_kernel_immediate;

template<>
struct kernel_traits<poisson_G_kernel_immediate>
{
	typedef build<location<space_3d> >::type test_input_t;
	typedef build<location<space_3d> >::type trial_input_t;
	typedef typename poisson_g_wall_immediate<space_3d>::type output_t;
	typedef space_3d::scalar_t result_t;
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	static unsigned const singularity_order = 1;
	static unsigned const singular_quadrature_order = 7;
};

class poisson_G_kernel_immediate :
	public kernel_base<poisson_G_kernel_immediate>,
	public reciprocal_distance_kernel<poisson_G_kernel_immediate>
{
public:
	typedef kernel_base<poisson_G_kernel> base_t;
	typedef base_t::test_input_t test_input_t;
	typedef base_t::trial_input_t trial_input_t;
	typedef base_t::scalar_t scalar_t;

	template <class test_input_t, class trial_input_t>
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return 1.0 / (x.get_x() - y.get_x()).norm() / (4.0 * M_PI);
	}
};


template <class scalar>
struct poisson_h_brick
{
	template <class wall>
	class brick : public wall
	{
	public:
		template <class test_input_t, class trial_input_t>
		brick(test_input_t const &test_input, trial_input_t const &trial_input) :
			wall(test_input, trial_input),
			m_poisson_h(-1.0 * wall::get_poisson_g())
		{
			auto r = wall::get_distance();
			auto rdn = wall::get_distance_vector().dot(trial_input.get_unit_normal()) / r;
			m_poisson_h *= rdn / r;
		}
		
		scalar const & get_poisson_h(void) const
		{
			return m_poisson_h;
		}
		
		scalar const & get_result(void) const
		{
			return m_poisson_h;
		}
		
	private:
		scalar m_poisson_h;
	};
};


template <class space>
struct poisson_h_wall : glue<
	poisson_h_brick<typename space::scalar_t>::template brick,
	typename poisson_g_wall<space>::type
> {};


// forward declaration
class poisson_H_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_H_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef typename poisson_h_wall<space_3d>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$ */
class poisson_H_kernel :
	public kernel_base<poisson_H_kernel>,
	public reciprocal_distance_kernel<poisson_H_kernel>
{
public:
	/** \brief the crtp base's type */
	typedef kernel_base<poisson_H_kernel> base_t;
	/** \brief type of the test input */
	typedef base_t::test_input_t test_input_t;
	/** \brief type of the trial input */
	typedef base_t::trial_input_t trial_input_t;
	/** \brief type of the scalar of the inputs */
	typedef base_t::scalar_t scalar_t;

	/** \brief evaluate kernel between test and trial positions
	* \param [in] x the test input
	* \param [in] y the trial input
	* \return the kernel value K(x,y)
	*/
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		// instantiate output
		typedef typename kernel_traits<poisson_H_kernel>::output_t output_t;
		output_t output(x, y);
		// select result from output
		return output.get_result();
	}
};

#endif // POISSON_KERNEL_HPP_INCLUDED


