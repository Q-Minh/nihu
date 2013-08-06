/**
 * \file poisson_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Poisson equation \f$\nabla^2 p = 0\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef POISSON_KERNEL_HPP_INCLUDED
#define POISSON_KERNEL_HPP_INCLUDED

#include <cmath> // log
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


/** \brief a brick representing a 3D Poisson kernel \f$1/4\pi r\f$
 * \tparam scalar the scalar of the coordinate space the distance is defined over
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
			m_poisson_g(-log(wall::get_distance()) / (2.0 * M_PI))
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


/** \brief a brick representing a 3D Poisson kernel \f$1/4\pi r\f$
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


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_g_brick into a wall
 * \tparam space the coordinate space the Poisson kernel is defined over
 */
template <class scalar>
struct poisson_2d_g_wall : build<
	distance_vector_brick<space<scalar, 2> >,
	distance_brick<scalar>,
	poisson_2d_g_brick<scalar>
> {};


/** \brief combination of ::distance_vector_brick, ::distance_brick and ::poisson_g_brick into a wall
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

/** \brief traits of the Poisson G kernel */
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

/** \brief 3D Poisson kernel \f$-\ln r/2\pi \f$ */
class poisson_2d_G_kernel :
	public kernel_base<poisson_2d_G_kernel>,
	public reciprocal_distance_kernel<poisson_2d_G_kernel>
{
public:
	using reciprocal_distance_kernel<poisson_2d_G_kernel>::estimate_complexity;
};


/** \brief 3D Poisson kernel \f$1/4\pi r\f$ */
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
			m_poisson_h(-1.0 * wall::get_poisson_g())
		{
			auto r = wall::get_distance();
			auto rdn = wall::get_distance_vector().dot(trial_input.get_unit_normal()) / r;
			m_poisson_h *= rdn / r;
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
struct poisson_3d_h_wall : glue<
	poisson_3d_h_brick<scalar>::template brick,
	typename poisson_3d_g_wall<scalar>::type
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
	static bool const is_symmetric = true;
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


typedef poisson_3d_G_kernel poisson_3d_SLP_kernel;
typedef poisson_2d_G_kernel poisson_2d_SLP_kernel;
typedef poisson_3d_H_kernel poisson_3d_DLP_kernel;


/** \brief analytical expression of the collocational singular integral over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class collocational_singular_integral_shortcut<
	poisson_3d_SLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename TrialField::nset_t,
			tria_0_shape_set
		>::value
	>::type
>
{
public:
	template <class result_t>
	static result_t &eval(
		result_t &result,
		poisson_3d_SLP_kernel const &,
		TestField const &test_field,
		TrialField const &trial_field)
	{
		auto const &C_old = trial_field.get_elem().get_coords();
		auto const &x0 = test_field.get_elem().get_center();
		static unsigned N = tria_1_elem::num_nodes;

		double r[N];	// radius lengths
		typename tria_1_elem::coords_t R, C;
		for (unsigned i = 0; i < N; ++i)
		{
			R.col(i) = C_old.col(i) - x0;
			r[i] = R.col(i).norm();
			R.col(i) /= r[i];
			C.col(i) = C_old.col(i) - C_old.col((i+1) % N);
			C.col(i) /= C.col(i).norm();
		}

		for (unsigned i = 0; i < N; ++i)
		{
			double theta = std::acos(R.col(i).dot(R.col((i+1) % N)));
			double alpha = std::acos(R.col(i).dot(C.col(i)));
			result(0,0) += r[i] * std::sin(alpha) * std::log(std::tan((alpha+theta)/2.0)/tan(alpha/2.0));
		}

		result(0,0) /= (4.0 * M_PI);

		return result;
	}
};


/** \brief analytical expression of the collocational singular integral over a constant triangle
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField>
class collocational_singular_integral_shortcut<
	poisson_3d_DLP_kernel, TestField, TrialField,
	typename std::enable_if<
		std::is_same<
			typename TrialField::nset_t,
			tria_0_shape_set
		>::value
	>::type
>
{
public:
	template <class result_t>
	constexpr static result_t &eval(
		result_t &result,
		poisson_3d_DLP_kernel const &,
		TestField const &,
		TrialField const &)
	{
		return result;
	}
};


#endif // POISSON_KERNEL_HPP_INCLUDED

