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


/** \brief a brick representing a distance vector \f${\bf r} = {\bf y} - {\bf x}\f$
 * \tparam space the coordinate space the distance is defined over
 */
template <class space>
struct distance_vector_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_distance_vector(trial_input.get_x()-test_input.get_x())
		{
		}
		
		/** \brief return distance vector
		 * \return distance vector
		 */
		typename space::location_t const &get_distance_vector(void) const
		{
			return m_distance_vector;
		}
		
	private:
		typename space::location_t m_distance_vector;
	};
};


/** \brief a brick representing a scalar distance \f$r = |{\bf r}|\f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct distance_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_distance(wall::get_distance_vector().norm())
		{
		}
		
		/** \brief return distance
		 * \return scalar distance
		 */
		scalar const & get_distance(void) const
		{
			return m_distance;
		}
		
	private:
		scalar m_distance;
	};
};


/** \brief a brick representing a poisson kernel \f$1/4\pi r\f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_g_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_poisson_g(1.0 / wall::get_distance() / (4.0 * M_PI))
		{
		}
		
		/** \brief return poisson g kernel
		 * \return poisson g kernel
		 */
		scalar const & get_poisson_g(void) const
		{
			return m_poisson_g;
		}
		
		/** \brief return poisson g kernel
		 * \return poisson g kernel
		 */
		scalar const & get_result(void) const
		{
			return m_poisson_g;
		}
		
	private:
		scalar m_poisson_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_g_brick into a wall
 * \tparam space the coordinate space the poisson kernel is defined over
 */
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
	typedef poisson_g_wall<space_3d>::type output_t;
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
	using reciprocal_distance_kernel<poisson_G_kernel>::estimate_complexity;
};


/*
class poisson_G_kernel_immediate;

template<>
struct kernel_traits<poisson_G_kernel_immediate>
{
	typedef build<location<space_3d> >::type test_input_t;
	typedef build<location<space_3d> >::type trial_input_t;
	typedef empty_wall output_t;
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

	template <class tsi_t, class tri_t>
	result_t operator()(tsi_t const &x, tri_t const &y) const
	{
		return 1.0 / (x.get_x() - y.get_x()).norm() / (4.0 * M_PI);
	}
};
*/


/** \brief a brick representing a poisson derivative kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_h_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel the kernel instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_t const &kernel) :
			wall(test_input, trial_input, kernel),
			m_poisson_h(-1.0 * wall::get_poisson_g())
		{
			auto r = wall::get_distance();
			auto rdn = wall::get_distance_vector().dot(trial_input.get_unit_normal()) / r;
			m_poisson_h *= rdn / r;
		}
		
		/** \brief return poisson h kernel
		 * \return poisson h kernel
		 */
		scalar const &get_poisson_h(void) const
		{
			return m_poisson_h;
		}
		
		/** \brief return poisson h kernel
		 * \return poisson h kernel
		 */
		scalar const & get_result(void) const
		{
			return get_poisson_h();
		}
		
	private:
		scalar m_poisson_h;
	};
};


/** \brief combination of poisson_g_wall and poisson_h_brick into a wall
 * \tparam space the coordinate space the poisson kernel is defined over
 */
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
	typedef poisson_h_wall<space_3d>::type output_t;
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
	using reciprocal_distance_kernel<poisson_H_kernel>::estimate_complexity;
};


typedef poisson_G_kernel poisson_SLP_kernel;
typedef poisson_H_kernel poisson_DLP_kernel; 

#endif // POISSON_KERNEL_HPP_INCLUDED


