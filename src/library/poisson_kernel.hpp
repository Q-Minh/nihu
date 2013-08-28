/**
 * \file poisson_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Poisson equation \f$ \nabla^2 p = 0 \f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef POISSON_KERNEL_HPP_INCLUDED
#define POISSON_KERNEL_HPP_INCLUDED

#include <cmath>
#include "../bem/global_definitions.hpp"
#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"
#include "../bem/interval_estimator.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"
#include "reciprocal_kernel_intervals.hpp"


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
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
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


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_2d_g_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_g_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	poisson_2d_g_brick<scalar>
> {};


// forward declaration
class poisson_2d_SLP_kernel;

/** \brief traits of the Poisson 2D G kernel */
template<>
struct kernel_traits<poisson_2d_SLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_2d_g_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 2D Poisson kernel \f$ -\ln r/2\pi \f$ */
class poisson_2d_SLP_kernel :
	public kernel_base<poisson_2d_SLP_kernel>
{
};



/** \brief a brick representing a 2D Poisson derivative kernel \f$ -1/2\pi r \cdot r'_{n_y}\f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_2d_h_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_h(-wall::get_rdny()/wall::get_distance() / (2.0 * M_PI))
		{
		}

		/** \brief return Poisson h kernel
		 * \return Poisson h kernel
		 */
		scalar const & get_poisson_h(void) const
		{
			return m_poisson_h;
		}

		/** \brief return Poisson h kernel
		 * \return Poisson h kernel
		 */
		scalar const & get_result(void) const
		{
			return m_poisson_h;
		}

	private:
		scalar m_poisson_h;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdny_brick and ::poisson_2d_h_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_h_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	poisson_2d_h_brick<scalar>
> {};


// forward declaration
class poisson_2d_DLP_kernel;

/** \brief traits of the Poisson 2D H kernel */
template<>
struct kernel_traits<poisson_2d_DLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_2d_h_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_2d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 2D Poisson kernel \f$ -1/2\pi r \cdot r'_{n_y} \f$ */
class poisson_2d_DLP_kernel :
	public kernel_base<poisson_2d_DLP_kernel>
{
};


/** \brief a brick representing a 2D Poisson derivative kernel \f$ -1/2\pi r \cdot r'_{n_x}\f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_2d_ht_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_ht(-wall::get_rdnx()/wall::get_distance() / (2.0 * M_PI))
		{
		}

		/** \brief return Poisson h kernel
		 * \return Poisson ht kernel
		 */
		scalar const & get_poisson_ht(void) const
		{
			return m_poisson_ht;
		}

		/** \brief return Poisson ht kernel
		 * \return Poisson ht kernel
		 */
		scalar const & get_result(void) const
		{
			return m_poisson_ht;
		}

	private:
		scalar m_poisson_ht;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick and ::poisson_2d_ht_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_ht_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	poisson_2d_ht_brick<scalar>
> {};


// forward declaration
class poisson_2d_DLPt_kernel;

/** \brief traits of the Poisson 2D Ht kernel */
template<>
struct kernel_traits<poisson_2d_DLPt_kernel>
{
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type test_input_t;
	/** \brief kernel test input type */
	typedef build<location<space_2d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_2d_ht_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_2d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 2D Poisson kernel \f$ -1/2\pi r \cdot r'_{n_x} \f$ */
class poisson_2d_DLPt_kernel :
	public kernel_base<poisson_2d_DLPt_kernel>
{
};



/** \brief a brick representing a 2D Poisson hypersingular kernel \f$ \dots \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_2d_hyper_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_hyper(
				1.0/(2.0 * M_PI)/wall::get_distance()/wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					2.0 * wall::get_rdny()*wall::get_rdnx()))
		{
		}

		/** \brief return Poisson hypersingular kernel
		 * \return Poisson hypersingular kernel
		 */
		result_t const & get_poisson_hyper(void) const
		{
			return m_poisson_hyper;
		}

		/** \brief return Poisson hypersingular kernel
		 * \return Poisson hypersingular kernel
		 */
		result_t const & get_result(void) const
		{
			return m_poisson_hyper;
		}

	private:
		result_t m_poisson_hyper;
	};
};


/** \brief combination of several bricks into a poisson_2d_hyper wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_hyper_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	rdny_brick<scalar>,
	rdnx_brick<scalar>,
	poisson_2d_hyper_brick<scalar>
> {};


// forward declaration
class poisson_2d_HSP_kernel;

/** \brief traits of the Poisson 2D Hypersingular kernel */
template<>
struct kernel_traits<poisson_2d_HSP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_2d>, normal_jacobian<space_2d>  >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_2d_hyper_wall<space_2d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 2D Poisson kernel \f$ \dots \f$ */
class poisson_2d_HSP_kernel :
	public kernel_base<poisson_2d_HSP_kernel>
{
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
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
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
class poisson_3d_SLP_kernel;

/** \brief traits of the Poisson 3D G kernel */
template<>
struct kernel_traits<poisson_3d_SLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_3d_g_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 1;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 3D Poisson kernel \f$ 1/4\pi r\f$ */
class poisson_3d_SLP_kernel :
	public kernel_base<poisson_3d_SLP_kernel>
{
};


/** \brief a brick representing a Poisson derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_y} \f$
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
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_h(-wall::get_poisson_g() * wall::get_rdny() / wall::get_distance())
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
	rdny_brick<scalar>,
	poisson_3d_g_brick<scalar>,
	poisson_3d_h_brick<scalar>
> {};


// forward declaration
class poisson_3d_DLP_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_3d_DLP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_3d_h_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef typename output_t::result_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 3D Poisson kernel \f$ -1/4\pi r^2 \cdot dr/dn \f$ */
class poisson_3d_DLP_kernel :
	public kernel_base<poisson_3d_DLP_kernel>
{
};


/** \brief a brick representing a Poisson derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_x} \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_3d_ht_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_ht(-wall::get_poisson_g() * wall::get_rdnx() / wall::get_distance())
		{
		}

		/** \brief return Poisson ht kernel
		 * \return Poisson ht kernel
		 */
		scalar const &get_poisson_ht(void) const
		{
			return m_poisson_ht;
		}

		/** \brief return Poisson ht kernel
		 * \return Poisson ht kernel
		 */
		scalar const & get_result(void) const
		{
			return get_poisson_ht();
		}

	private:
		scalar m_poisson_ht;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick, poisson_3d_g_brick and poisson_3d_ht_brick into a wall
 * \tparam scalar the scalar type of the poisson ht kernel result
 */
template <class scalar>
struct poisson_3d_ht_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	poisson_3d_g_brick<scalar>,
	poisson_3d_ht_brick<scalar>
> {};


// forward declaration
class poisson_3d_DLPt_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_3d_DLPt_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_3d_ht_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = false;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 2;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 3D Poisson derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_x} \f$ */
class poisson_3d_DLPt_kernel :
	public kernel_base<poisson_3d_DLPt_kernel>
{
};


/** \brief a brick representing a Poisson double derivative kernel \f$ 1/4\pi r^3 \cdot \left( n_x n_y + 3 r'_{n_x} r'_{n_y} \right) \f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar>
struct poisson_3d_hyper_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef scalar result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel type
		 * \param [in] test_input the test input
		 * \param [in] trial_input the trial input
		 * \param [in] kernel_data the kernel data instance
		 */
		template <class test_input_t, class trial_input_t, class kernel_data_t>
		brick(
			test_input_t const &test_input,
			trial_input_t const &trial_input,
			kernel_data_t const &kernel_data) :
			wall(test_input, trial_input, kernel_data),
			m_poisson_hyper(
				wall::get_poisson_g() / wall::get_distance() / wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					3.0 * wall::get_rdnx() * wall::get_rdny()
				))
		{
		}

		/** \brief return Poisson ht kernel
		 * \return Poisson ht kernel
		 */
		scalar const &get_poisson_hyper(void) const
		{
			return m_poisson_hyper;
		}

		/** \brief return Poisson ht kernel
		 * \return Poisson ht kernel
		 */
		scalar const & get_result(void) const
		{
			return get_poisson_hyper();
		}

	private:
		scalar m_poisson_hyper;
	};
};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick, poisson_3d_g_brick and poisson_3d_ht_brick into a wall
 * \tparam scalar the scalar type of the poisson ht kernel result
 */
template <class scalar>
struct poisson_3d_hyper_wall : build<
	distance_vector_brick<space<scalar, 3> >,
	distance_brick<scalar>,
	rdnx_brick<scalar>,
	rdny_brick<scalar>,
	poisson_3d_g_brick<scalar>,
	poisson_3d_hyper_brick<scalar>
> {};


// forward declaration
class poisson_3d_HSP_kernel;

/** \brief traits of the Poisson H kernel */
template<>
struct kernel_traits<poisson_3d_HSP_kernel>
{
	/** \brief kernel test input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type test_input_t;
	/** \brief kernel trial input type */
	typedef build<location<space_3d>, normal_jacobian<space_3d> >::type trial_input_t;
	/** \brief the data type */
	typedef collect<empty_data> data_t;
	/** \brief the kernel output type */
	typedef poisson_3d_hyper_wall<space_3d::scalar_t>::type output_t;
	/** \brief kernel result type */
	typedef space_3d::scalar_t result_t;
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief kernel singularity order ( r^(-order) ) */
	static unsigned const singularity_order = 3;
	/** \brief quadrature order used to generate Duffy singular quadratures */
	static unsigned const singular_quadrature_order = 7;
	/** \brief the kernel complexity estimator class */
	typedef interval_estimator<
		typename reciprocal_distance_kernel_interval<singularity_order, GLOBAL_ACCURACY>::type
	> complexity_estimator_t;
};

/** \brief 3D Poisson derivative kernel \f$ 1/4\pi r^3 \cdot \left( n_x n_y + 3 r'_{n_x} r'_{n_y} \right) \f$ */
class poisson_3d_HSP_kernel :
	public kernel_base<poisson_3d_HSP_kernel>
{
};

#endif // POISSON_KERNEL_HPP_INCLUDED

