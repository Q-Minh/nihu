// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file laplace_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Laplace equation \f$ \nabla^2 p = 0 \f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef LAPLACE_KERNEL_HPP_INCLUDED
#define LAPLACE_KERNEL_HPP_INCLUDED

#include <cmath>
#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../library/potential_kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"

/** \brief a brick representing a Laplace kernel
 * \tparam Space the coordinate space the kernel is defined over
 */
template <class Space, class Layer>
struct laplace_brick;

/** \brief combination of several bricks into a Laplace wall
 * \tparam Space the coordinate space the Laplace kernel is defined over
 */
template <class Space, class Layer>
struct laplace_wall;

/** \brief kernel of the Laplace equation
 * \tparam Space the coordinate space the kernel is defined over
 */
template <class Space, class Layer>
class laplace_kernel : public kernel_base<laplace_kernel<Space, Layer> >
{
};


/** \brief specialisation of ::laplace_brick for the 2D SLP case */
template <class Scalar>
struct laplace_brick<space_2d<Scalar>, potential::SLP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_g(-std::log(wall::get_distance()) / (2. * M_PI))
		{
		}

		/** \brief return Laplace g kernel
		 * \return Laplace g kernel
		 */
		result_t const & get_laplace_g(void) const
		{
			return m_laplace_g;
		}

		/** \brief return Laplace g kernel
		 * \return Laplace g kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_g;
		}

	private:
		result_t m_laplace_g;
	};
};

/** \brief specialisation of ::laplace_brick for the 3D SLP case */
template <class Scalar>
struct laplace_brick<space_3d<Scalar>, potential::SLP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_g(1. / wall::get_distance() / (4. * M_PI))
		{
		}

		/** \brief return Laplace g kernel
		 * \return Laplace g kernel
		 */
		result_t const & get_laplace_g(void) const
		{
			return m_laplace_g;
		}

		/** \brief return Laplace g kernel
		 * \return Laplace g kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_g;
		}

	private:
		result_t m_laplace_g;
	};
};

/** \brief specialisation of ::laplace_wall for the SLP case */
template <class Space>
struct laplace_wall<Space, potential::SLP> : build<
	distance_vector_brick<Space>,
	distance_brick<typename Space::scalar_t>,
	laplace_brick<Space, potential::SLP>
> {};

namespace kernel_traits_ns
{
	template <class Space, class Layer>
	struct space<laplace_kernel<Space, Layer> > : Space {};

	template <class Space>
	struct test_input<laplace_kernel<Space, potential::SLP> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::SLP> > : build<location<Space> > {};

	template <class Space, class Layer>
	struct data<laplace_kernel<Space, Layer> >
	{
		typedef collect<empty_data> type;
	};

	template <class Space, class Layer>
	struct output<laplace_kernel<Space, Layer> > : laplace_wall<Space, Layer> {};

	template <class Space, class Layer>
	struct quadrature_family<laplace_kernel<Space, Layer> > : gauss_family_tag {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::SLP> > : asymptotic::log<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::SLP> > : asymptotic::inverse<1> {};

	template <class Space, class Layer>
	struct result_dimension<laplace_kernel<Space, Layer> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::SLP> > : std::true_type {};

	template <class Space, class Layer>
	struct is_singular<laplace_kernel<Space, Layer> > : std::true_type {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::SLP> > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::SLP> > : asymptotic::inverse<1> {};

	/** \brief the singular quadrature order of the laplace SLP kernel
	 * \todo check if the same value can be used for the 2D and 3D case
	 */
	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::SLP> > : std::integral_constant<unsigned, 7> {};

	template <class Space, class Layer>
	struct singular_core<laplace_kernel<Space, Layer> > {
		typedef  laplace_kernel<Space, Layer>  type;
	};
}

/** \brief specialisation of ::laplace_brick for the 2D DLP case */
template <class Scalar>
struct laplace_brick<space_2d<Scalar>, potential::DLP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_h(-wall::get_rdny()/wall::get_distance() / (2. * M_PI))
		{
		}

		/** \brief return Laplace h kernel
		 * \return Laplace h kernel
		 */
		result_t const & get_laplace_h(void) const
		{
			return m_laplace_h;
		}

		/** \brief return Laplace h kernel
		 * \return Laplace h kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_h;
		}

	private:
		result_t m_laplace_h;
	};
};

/** \brief specialisation of ::laplace_brick for the 3D DLP case */
template <class Scalar>
struct laplace_brick<space_3d<Scalar>, potential::DLP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_h(-wall::get_laplace_g() * wall::get_rdny() / wall::get_distance())
		{
		}

		/** \brief return Laplace h kernel
		 * \return Laplace h kernel
		 */
		result_t const &get_laplace_h(void) const
		{
			return m_laplace_h;
		}

		/** \brief return Laplace h kernel
		 * \return Laplace h kernel
		 */
		result_t const & get_result(void) const
		{
			return get_laplace_h();
		}

	private:
		result_t m_laplace_h;
	};
};

/** \brief specialisation of ::laplace_wall for the 2D DLP case */
template <class Scalar>
struct laplace_wall<space_2d<Scalar>, potential::DLP> : build<
	distance_vector_brick<space_2d<Scalar>>,
	distance_brick<Scalar>,
	rdny_brick<Scalar>,
	laplace_brick<space_2d<Scalar>, potential::DLP>
> {};

/** \brief specialisation of ::laplace_wall for the 3D DLP case */
template <class Scalar>
struct laplace_wall<space_3d<Scalar>, potential::DLP> : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	rdny_brick<Scalar>,
	laplace_brick<space_3d<Scalar>, potential::SLP>,
	laplace_brick<space_3d<Scalar>, potential::DLP>
> {};

namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::DLP> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::DLP> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::DLP> > : std::false_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::DLP> > : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::DLP> > : asymptotic::inverse<2> {};

	/** \brief kernel singularity type
	 * \todo check this!
	 */
	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::DLP> > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::DLP> > : asymptotic::inverse<1> {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::DLP> > : std::integral_constant<unsigned, 7> {};
}

/** \brief specialisation of ::laplace_brick for the 2D DLPt case */
template <class Scalar>
struct laplace_brick<space_2d<Scalar>, potential::DLPt>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_ht(-wall::get_rdnx()/wall::get_distance() / (2. * M_PI))
		{
		}

		/** \brief return Laplace h kernel
		 * \return Laplace ht kernel
		 */
		result_t const & get_laplace_ht(void) const
		{
			return m_laplace_ht;
		}

		/** \brief return Laplace ht kernel
		 * \return Laplace ht kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_ht;
		}

	private:
		result_t m_laplace_ht;
	};
};

/** \brief specialisation of ::laplace_brick for the 3D DLPt case */
template <class Scalar>
struct laplace_brick<space_3d<Scalar>, potential::DLPt>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_ht(-wall::get_laplace_g() * wall::get_rdnx() / wall::get_distance())
		{
		}

		/** \brief return Laplace ht kernel
		 * \return Laplace ht kernel
		 */
		result_t const &get_laplace_ht(void) const
		{
			return m_laplace_ht;
		}

		/** \brief return Laplace ht kernel
		 * \return Laplace ht kernel
		 */
		result_t const & get_result(void) const
		{
			return get_laplace_ht();
		}

	private:
		result_t m_laplace_ht;
	};
};

/** \brief specialsation of ::laplace_wall for the 2D DLPt case */
template <class Scalar>
struct laplace_wall<space_2d<Scalar>, potential::DLPt> : build<
	distance_vector_brick<space_2d<Scalar> >,
	distance_brick<Scalar>,
	rdnx_brick<Scalar>,
	laplace_brick<space_2d<Scalar>, potential::DLPt>
> {};

/** \brief specialsation of ::laplace_wall for the 3D DLPt case */
template <class Scalar>
struct laplace_wall<space_3d<Scalar>, potential::DLPt> : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	rdnx_brick<Scalar>,
	laplace_brick<space_3d<Scalar>, potential::SLP>,
	laplace_brick<space_3d<Scalar>, potential::DLPt>
> {};

namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::DLPt> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::DLPt> > : build<location<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::DLPt> > : std::false_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::DLPt> > : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::DLPt> > : asymptotic::inverse<2> {};

	/** \brief the singularity type
	 * \todo check this just like the plain DLP kernel
	 */
	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::DLPt> > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::DLPt> > : asymptotic::inverse<1> {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::DLPt> > : std::integral_constant<unsigned, 7> {};
}

/** \brief specialisation of ::laplace_brick for the 2D HSP case */
template <class Scalar>
struct laplace_brick<space_2d<Scalar>, potential::HSP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_hyper(
				1./(2. * M_PI)/wall::get_distance()/wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					2. * wall::get_rdny()*wall::get_rdnx()))
		{
		}

		/** \brief return Laplace hypersingular kernel
		 * \return Laplace hypersingular kernel
		 */
		result_t const & get_laplace_hyper(void) const
		{
			return m_laplace_hyper;
		}

		/** \brief return Laplace hypersingular kernel
		 * \return Laplace hypersingular kernel
		 */
		result_t const & get_result(void) const
		{
			return m_laplace_hyper;
		}

	private:
		result_t m_laplace_hyper;
	};
};

/** \brief specialisation of ::laplace_brick for the 3D HSP case */
template <class Scalar>
struct laplace_brick<space_3d<Scalar>, potential::HSP>
{
	/** \brief the brick template
	 * \tparam wall the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef Scalar result_t;

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
			m_laplace_hyper(
				wall::get_laplace_g() / wall::get_distance() / wall::get_distance() * (
					test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					3. * wall::get_rdnx() * wall::get_rdny()
				))
		{
		}

		/** \brief return Laplace hypersingular kernel
		 * \return Laplace hypersingular kernel
		 */
		result_t const &get_laplace_hyper(void) const
		{
			return m_laplace_hyper;
		}

		/** \brief return Laplace hypersingular kernel
		 * \return Laplace hypersingular kernel
		 */
		result_t const & get_result(void) const
		{
			return get_laplace_hyper();
		}

	private:
		result_t m_laplace_hyper;
	};
};

/** \brief specialisation of ::laplace_wall for the 2D HSP case */
template <class Scalar>
struct laplace_wall<space_2d<Scalar>, potential::HSP> : build<
	distance_vector_brick<space_2d<Scalar> >,
	distance_brick<Scalar>,
	rdny_brick<Scalar>,
	rdnx_brick<Scalar>,
	laplace_brick<space_2d<Scalar>, potential::HSP>
> {};

/** \brief specialisation of ::laplace_wall for the 3D HSP case */
template <class Scalar>
struct laplace_wall<space_3d<Scalar>, potential::HSP> : build<
	distance_vector_brick<space_3d<Scalar>>,
	distance_brick<Scalar>,
	rdnx_brick<Scalar>,
	rdny_brick<Scalar>,
	laplace_brick<space_3d<Scalar>, potential::SLP>,
	laplace_brick<space_3d<Scalar>, potential::HSP>
> {};

namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::HSP> > : build<location<Space>, normal_jacobian<Space > > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::HSP> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::HSP> > : std::true_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::HSP> > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::HSP> > : asymptotic::inverse<3> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::HSP> > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::HSP> > : asymptotic::inverse<3> {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::HSP> > : std::integral_constant<unsigned, 7> {};
}

/** \brief shorthand for the 2d laplace SLP kernel */
typedef laplace_kernel<space_2d<>, potential::SLP> laplace_2d_SLP_kernel;
/** \brief shorthand for the 3d laplace SLP kernel */
typedef laplace_kernel<space_3d<>, potential::SLP> laplace_3d_SLP_kernel;
/** \brief shorthand for the 2d laplace DLP kernel */
typedef laplace_kernel<space_2d<>, potential::DLP> laplace_2d_DLP_kernel;
/** \brief shorthand for the 3d laplace DLP kernel */
typedef laplace_kernel<space_3d<>, potential::DLP> laplace_3d_DLP_kernel;
/** \brief shorthand for the 2d laplace DLPt kernel */
typedef laplace_kernel<space_2d<>, potential::DLPt> laplace_2d_DLPt_kernel;
/** \brief shorthand for the 3d laplace DLPt kernel */
typedef laplace_kernel<space_3d<>, potential::DLPt> laplace_3d_DLPt_kernel;
/** \brief shorthand for the 2d laplace HSP kernel */
typedef laplace_kernel<space_2d<>, potential::HSP> laplace_2d_HSP_kernel;
/** \brief shorthand for the 3d laplace HSP kernel */
typedef laplace_kernel<space_3d<>, potential::HSP> laplace_3d_HSP_kernel;


#include "guiggiani_1992.hpp"

/** \brief specialisation of class ::polar_laurent_coeffs for the ::laplace_3d_HSP_kernel */
template <class Scalar>
class polar_laurent_coeffs<laplace_kernel<space_3d<Scalar>, potential::HSP> >
{
public:
	template <class guiggiani>
	static void eval(guiggiani &obj)
	{
		auto g1vec = obj.get_rvec_series(_1()) * (
			obj.get_rvec_series(_2()).dot(obj.get_Jvec_series(_0()))
			+ obj.get_rvec_series(_1()).dot(obj.get_Jvec_series(_1()))
			);

		auto b0vec = -obj.get_Jvec_series(_0());
		auto b1vec = 3. * g1vec - obj.get_Jvec_series(_1());

		auto a0 = b0vec.dot(obj.get_n0()) * obj.get_shape_series(_0());
		auto a1 = b1vec.dot(obj.get_n0()) * obj.get_shape_series(_0())
			+ b0vec.dot(obj.get_n0()) * obj.get_shape_series(_1());

		auto Sm2 = -3. * obj.get_rvec_series(_1()).dot(obj.get_rvec_series(_2()));

		obj.set_laurent_coeff(_m1(), -(Sm2 * a0 + a1) / (4. * M_PI));
		obj.set_laurent_coeff(_m2(), -a0 / (4. * M_PI));
	}
};

#endif // LAPLACE_KERNEL_HPP_INCLUDED

