/**
 * \file poisson_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Poisson equation \f$ \nabla^2 p = 0 \f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef POISSON_KERNEL_HPP_INCLUDED
#define POISSON_KERNEL_HPP_INCLUDED

#include <cmath>
#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

#include "location_normal.hpp"
#include "reciprocal_distance_kernel.hpp"

#include "basic_bricks.hpp"

/** \brief a brick representing a 2D Poisson kernel \f$ -\log r /2\pi \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_2d_g_brick
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
			m_poisson_g(-std::log(wall::get_distance()) / (2.0 * M_PI))
		{
		}

		/** \brief return Poisson g kernel
		 * \return Poisson g kernel
		 */
		scalar const & get_poisson_g(void) const
		{
			return m_poisson_g;
		}

		/** \brief return Poisson g kernel
		 * \return Poisson g kernel
		 */
		scalar const & get_result(void) const
		{
			return m_poisson_g;
		}

	private:
		scalar m_poisson_g;
	};
};


/** \brief a brick representing a 3D Poisson kernel \f$ 1/4\pi r \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_3d_g_brick
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

		/** \brief return Poisson g kernel
		 * \return Poisson g kernel
		 */
		scalar const & get_poisson_g(void) const
		{
			return m_poisson_g;
		}

		/** \brief return Poisson g kernel
		 * \return Poisson g kernel
		 */
		scalar const & get_result(void) const
		{
			return m_poisson_g;
		}

	private:
		scalar m_poisson_g;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_2d_g_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_g_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	poisson_2d_g_brick<scalar>
> {};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_3d_g_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_3d_g_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	poisson_3d_g_brick<scalar>
> {};


// forward declaration
class poisson_2d_G_kernel;

// forward declaration
class poisson_3d_G_kernel;

/** \brief traits of the Poisson 2D G kernel */
template<>
struct kernel_traits<poisson_2d_G_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef poisson_2d_g_wall<double>::type output_t;
	/** \brief kernel result type */
	typedef double result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief traits of the Poisson 3D G kernel */
template<>
struct kernel_traits<poisson_3d_G_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef poisson_3d_g_wall<double>::type output_t;
	/** \brief kernel result type */
	typedef double result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 2D Poisson kernel \f$ -\ln r/2\pi \f$ */
class poisson_2d_G_kernel :
	public kernel_base<poisson_2d_G_kernel>,
	public reciprocal_distance_kernel<poisson_2d_G_kernel>
{
public:
	using reciprocal_distance_kernel<poisson_2d_G_kernel>::estimate_complexity;
};


/** \brief 3D Poisson kernel \f$ 1/4\pi r\f$ */
class poisson_3d_G_kernel :
	public kernel_base<poisson_3d_G_kernel>,
	public reciprocal_distance_kernel<poisson_3d_G_kernel>
{
public:
	using reciprocal_distance_kernel<poisson_3d_G_kernel>::estimate_complexity;
};


/** \brief a brick representing a Poisson derivative kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_3d_h_brick
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
			m_poisson_h(-1.0 * wall::get_poisson_g() * wall::get_rdny() / wall::get_distance())
		{
		}

		/** \brief return Poisson h kernel
		 * \return Poisson h kernel
		 */
		scalar const &get_poisson_h(void) const
		{
			return m_poisson_h;
		}

		/** \brief return Poisson h kernel
		 * \return Poisson h kernel
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
template <class scalar>
struct poisson_3d_h_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	poisson_3d_g_brick<scalar>,
	rdny_brick<scalar>,
	poisson_3d_h_brick<scalar>
> {};


// forward declaration
class poisson_3d_H_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_3d_H_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the kernel output type */
	typedef poisson_3d_h_wall<double>::type output_t;
	/** \brief kernel result type */
	typedef double result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
};

/** \brief 3D Poisson kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$ */
class poisson_3d_H_kernel :
	public kernel_base<poisson_3d_H_kernel>,
	public reciprocal_distance_kernel<poisson_3d_H_kernel>
{
public:
	using reciprocal_distance_kernel<poisson_3d_H_kernel>::estimate_complexity;
};


typedef poisson_2d_G_kernel poisson_2d_SLP_kernel;
// typedef poisson_2d_H_kernel poisson_2d_DLP_kernel;

typedef poisson_3d_G_kernel poisson_3d_SLP_kernel;
typedef poisson_3d_H_kernel poisson_3d_DLP_kernel;

#endif // POISSON_KERNEL_HPP_INCLUDED

