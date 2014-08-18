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
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"


/** \brief a brick representing a Laplace SLP kernel
 * \tparam Space the coordinate space ther kernel is defined over
 */
template <class Space>
struct laplace_SLP_brick;

template <class Scalar>
struct laplace_SLP_brick<space_2d<Scalar> >
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

template <class Scalar>
struct laplace_SLP_brick<space_3d<Scalar> >
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


/** \brief combination of several bricks into a Laplace SLP wall
 * \tparam Space the coordinate space the Laplace kernel is defined over
 */
template <class Space>
struct laplace_SLP_wall : build<
	distance_vector_brick<Space>,
	distance_brick<typename Space::scalar_t>,
	laplace_SLP_brick<Space>
> {};


/** \brief Single layer potential kernel of the Laplace equation
 * \tparam Space the coordinate space the kernel is defined over
 */
template <class Space>
class laplace_SLP_kernel;

namespace kernel_traits_ns
{
	template <class Space> 
	struct space<laplace_SLP_kernel<Space> > : Space {};

	template <class Space> 
	struct test_input<laplace_SLP_kernel<Space> > : build<location<Space> > {};

	template <class Space> 
	struct trial_input<laplace_SLP_kernel<Space> > : build<location<Space> > {};

	template <class Space> 
	struct data<laplace_SLP_kernel<Space> >
	{
		typedef collect<empty_data> type;
	};

	template <class Space> 
	struct output<laplace_SLP_kernel<Space> > : laplace_SLP_wall<Space> {};

	template <class Space> 
	struct quadrature_family<laplace_SLP_kernel<Space> > : gauss_family_tag {};

	template <class Scalar> 
	struct far_field_behaviour<laplace_SLP_kernel<space_2d<Scalar> > > : asymptotic::log<1> {};

	template <class Scalar> 
	struct far_field_behaviour<laplace_SLP_kernel<space_3d<Scalar> > > : asymptotic::inverse<1> {};

	template <class Space> 
	struct result_dimension<laplace_SLP_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space> 
	struct is_symmetric<laplace_SLP_kernel<Space> > : std::true_type {};

	template <class Space> 
	struct is_singular<laplace_SLP_kernel<Space> > : std::true_type {};

	template <class Scalar>
	struct singularity_type<laplace_SLP_kernel<space_2d<Scalar> > > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_SLP_kernel<space_3d<Scalar> > > : asymptotic::inverse<1> {};

	/** \brief the singular quadrature order of the laplace SLP kernel
	 * \todo check if the same value can be used for the 2D and 3D case
	 */
	template <class Space>
	struct singular_quadrature_order<laplace_SLP_kernel<Space> > : std::integral_constant<unsigned, 7> {};

	template <class Space>
	struct singular_core<laplace_SLP_kernel<Space> > {
		typedef  laplace_SLP_kernel<Space>  type;
	};
}


template <class Space>
class laplace_SLP_kernel :
	public kernel_base<laplace_SLP_kernel<Space> >
{
};





/** \brief a brick representing a Laplace normal derivative kernel \f$ -1/2\pi r \cdot r'_{n_y}\f$
 * \tparam Space the coordinate space the distance is defined over
 */
template <class Space>
struct laplace_DLP_brick;

template <class Scalar>
struct laplace_DLP_brick<space_2d<Scalar> >
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

/** \brief a brick representing a 3d Laplace DLP kernel \f$ -1/4\pi r^2 \cdot r'_{n_y} \f$
 * \tparam Scalar the scalar of the coordinate space the distance is defined over
 */
template <class Scalar>
struct laplace_DLP_brick<space_3d<Scalar> >
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


template <class Space>
struct laplace_DLP_wall;

template <class Space>
struct laplace_DLP_wall : build<
	distance_vector_brick<Space>,
	distance_brick<typename Space::scalar_t>,
	rdny_brick<typename Space::scalar_t>,
	laplace_DLP_brick<Space>
> {};

template <class Scalar>
struct laplace_DLP_wall<space_3d<Scalar> > : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	rdny_brick<Scalar>,
	laplace_SLP_brick<space_3d<Scalar> >,
	laplace_DLP_brick<space_3d<Scalar> >
> {};



// forward declaration
template <class Space>
class laplace_DLP_kernel;

/** \brief traits of the Laplace DLP kernel */
namespace kernel_traits_ns
{
	template <class Space> 
	struct space<laplace_DLP_kernel<Space> > : Space {};

	template <class Space>
	struct test_input<laplace_DLP_kernel<Space> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_DLP_kernel<Space> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct data<laplace_DLP_kernel<Space> >
	{
		typedef collect<empty_data> type;
	};

	template <class Space>
	struct output<laplace_DLP_kernel<Space> > : laplace_DLP_wall<Space> {};

	template <class Space>
	struct result_dimension<laplace_DLP_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct quadrature_family<laplace_DLP_kernel<Space> > : gauss_family_tag {};

	template <class Space>
	struct is_symmetric<laplace_DLP_kernel<Space> > : std::false_type {};

	template <class Space>
	struct is_singular<laplace_DLP_kernel<Space> > : std::true_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_DLP_kernel<space_2d<Scalar> > > : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_DLP_kernel<space_3d<Scalar> > > : asymptotic::inverse<2> {};

	/** \brief kernel singularity type
	 * \todo check this!
	 */
	template <class Scalar>
	struct singularity_type<laplace_DLP_kernel<space_2d<Scalar> > > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_DLP_kernel<space_3d<Scalar> > > : asymptotic::inverse<1> {};

	template <class Space>
	struct singular_quadrature_order<laplace_DLP_kernel<Space> > : std::integral_constant<unsigned, 7> {};

	template <class Space>
	struct singular_core<laplace_DLP_kernel<Space> > {
		typedef  laplace_DLP_kernel<Space>  type;
	};
}


/** \brief Double layer potential kernel of the Laplace equation */
template <class Space>
class laplace_DLP_kernel :
	public kernel_base<laplace_DLP_kernel<Space> >
{
};




/** \brief a brick representing a Laplace transpose DLP kernel
 * \tparam Space the coordinate space the distance is defined over
 */
template <class Space>
struct laplace_DLPt_brick;


template <class Scalar>
struct laplace_DLPt_brick<space_2d<Scalar> >
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

/** \brief a brick representing a Laplace derivative kernel \f$ -1/4\pi r^2 \cdot r'_{n_x} \f$
 * \tparam Scalar the scalar of the coordinate space the distance is defined over
 */
template <class Scalar>
struct laplace_DLPt_brick<space_3d<Scalar> >
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


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick and ::laplace_2d_DLPt_brick into a wall
 * \tparam Space the coordinate space the Laplace kernel is defined over
 */
template <class Space>
struct laplace_DLPt_wall : build<
	distance_vector_brick<Space>,
	distance_brick<typename Space::scalar_t>,
	rdnx_brick<typename Space::scalar_t>,
	laplace_DLPt_brick<Space>
> {};


/** \brief combination of ::distance_vector_brick, ::distance_brick, ::rdnx_brick, laplace_3d_SLP_brick and laplace_3d_DLPt_brick into a wall
 * \tparam Scalar the scalar type of the Laplace ht kernel result
 */
template <class Scalar>
struct laplace_DLPt_wall<space_3d<Scalar> > : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	rdnx_brick<Scalar>,
	laplace_SLP_brick<space_3d<Scalar> >,
	laplace_DLPt_brick<space_3d<Scalar> >
> {};


// forward declaration
template <class Space>
class laplace_DLPt_kernel;

namespace kernel_traits_ns
{
	template <class Space> 
	struct space<laplace_DLPt_kernel<Space> > : Space {};

	template <class Space>
	struct test_input<laplace_DLPt_kernel<Space> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct trial_input<laplace_DLPt_kernel<Space> > : build<location<Space> > {};

	template <class Space>
	struct data<laplace_DLPt_kernel<Space> > {
		typedef collect<empty_data> type;
	};

	template <class Space>
	struct output<laplace_DLPt_kernel<Space> > : laplace_DLPt_wall<Space> {};

	template <class Space>
	struct result_dimension<laplace_DLPt_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct quadrature_family<laplace_DLPt_kernel<Space> > : gauss_family_tag {};

	template <class Space>
	struct is_symmetric<laplace_DLPt_kernel<Space> > : std::false_type {};

	template <class Space>
	struct is_singular<laplace_DLPt_kernel<Space> > : std::true_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_DLPt_kernel<space_2d<Scalar> > > : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_DLPt_kernel<space_3d<Scalar> > > : asymptotic::inverse<2> {};

	/** \brief the singularity type
	 * \todo check this just like the plain DLP kernel
	 */
	template <class Scalar>
	struct singularity_type<laplace_DLPt_kernel<space_2d<Scalar> > > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_DLPt_kernel<space_3d<Scalar> > > : asymptotic::inverse<1> {};

	template <class Space>
	struct singular_quadrature_order<laplace_DLPt_kernel<Space> > : std::integral_constant<unsigned, 7> {};

	template <class Space>
	struct singular_core<laplace_DLPt_kernel<Space> > {
		typedef  laplace_DLPt_kernel<Space>  type;
	};
}


/** \brief Transposed DLP kernel of the Laplace equation */
template <class Space>
class laplace_DLPt_kernel :
	public kernel_base<laplace_DLPt_kernel<Space> >
{
};




template <class Space>
struct laplace_HSP_brick;

/** \brief a brick representing a 2D Laplace hypersingular kernel
 * \tparam Scalar scalar type of the coordinate space the kernel is defined over
 */
template <class Scalar>
struct laplace_HSP_brick<space_2d<Scalar> >
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


/** \brief a brick representing a Laplace double derivative kernel \f$ 1/4\pi r^3 \cdot \left( n_x n_y + 3 r'_{n_x} r'_{n_y} \right) \f$
 * \tparam Scalar the scalar of the coordinate space the distance is defined over
 */
template <class Scalar>
struct laplace_HSP_brick<space_3d<Scalar> >
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


/** \brief combination of several bricks into a laplace hypersingular wall
 * \tparam Scalar the coordinate space the Laplace kernel is defined over
 */
template <class Space>
struct laplace_HSP_wall;

template <class Scalar>
struct laplace_HSP_wall<space_2d<Scalar> > : build<
	distance_vector_brick<space_2d<Scalar> >,
	distance_brick<Scalar>,
	rdny_brick<Scalar>,
	rdnx_brick<Scalar>,
	laplace_HSP_brick<space_2d<Scalar> >
> {};

template <class Scalar>
struct laplace_HSP_wall<space_3d<Scalar> > : build<
	distance_vector_brick<space_3d<Scalar>>,
	distance_brick<Scalar>,
	rdnx_brick<Scalar>,
	rdny_brick<Scalar>,
	laplace_SLP_brick<space_3d<Scalar> >,
	laplace_HSP_brick<space_3d<Scalar> >
> {};


// forward declaration
template <class Space>
class laplace_HSP_kernel;

namespace kernel_traits_ns
{
	template <class Space> 
	struct space<laplace_HSP_kernel<Space> > : Space {};

	template <class Space>
	struct test_input<laplace_HSP_kernel<Space> > : build<location<Space>, normal_jacobian<Space > > {};

	template <class Space>
	struct trial_input<laplace_HSP_kernel<Space> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct data<laplace_HSP_kernel<Space> > {
		typedef collect<empty_data> type;
	};

	template <class Space>
	struct output<laplace_HSP_kernel<Space> > : laplace_HSP_wall<Space> {};

	template <class Space>
	struct result_dimension<laplace_HSP_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct quadrature_family<laplace_HSP_kernel<Space> > : gauss_family_tag {};

	template <class Space>
	struct is_symmetric<laplace_HSP_kernel<Space> > : std::true_type {};

	template <class Space>
	struct is_singular<laplace_HSP_kernel<Space> > : std::true_type {};

	template <class Scalar>
	struct far_field_behaviour<laplace_HSP_kernel<space_2d<Scalar> > > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_HSP_kernel<space_3d<Scalar> > > : asymptotic::inverse<3> {};

	template <class Scalar>
	struct singularity_type<laplace_HSP_kernel<space_2d<Scalar> > > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<laplace_HSP_kernel<space_3d<Scalar> > > : asymptotic::inverse<3> {};

	template <class Space>
	struct singular_quadrature_order<laplace_HSP_kernel<Space> > : std::integral_constant<unsigned, 7> {};

	template <class Space>
	struct singular_core<laplace_HSP_kernel<Space> > {
		typedef  laplace_HSP_kernel<Space>  type;
	};
}


/** \brief Hypersingular kernel of the Laplace equation */
template <class Space>
class laplace_HSP_kernel :
	public kernel_base<laplace_HSP_kernel<Space> >
{
};


/** \brief shorthand for the 2d laplace SLP kernel */
typedef laplace_SLP_kernel<space_2d<> > laplace_2d_SLP_kernel;
/** \brief shorthand for the 3d laplace SLP kernel */
typedef laplace_SLP_kernel<space_3d<> > laplace_3d_SLP_kernel;
/** \brief shorthand for the 2d laplace DLP kernel */
typedef laplace_DLP_kernel<space_2d<> > laplace_2d_DLP_kernel;
/** \brief shorthand for the 3d laplace DLP kernel */
typedef laplace_DLP_kernel<space_3d<> > laplace_3d_DLP_kernel;
/** \brief shorthand for the 2d laplace DLPt kernel */
typedef laplace_DLPt_kernel<space_2d<> > laplace_2d_DLPt_kernel;
/** \brief shorthand for the 3d laplace DLPt kernel */
typedef laplace_DLPt_kernel<space_3d<> > laplace_3d_DLPt_kernel;
/** \brief shorthand for the 2d laplace HSP kernel */
typedef laplace_HSP_kernel<space_2d<> > laplace_2d_HSP_kernel;
/** \brief shorthand for the 3d laplace HSP kernel */
typedef laplace_HSP_kernel<space_3d<> > laplace_3d_HSP_kernel;



#include "guiggiani_1992.hpp"

/** \brief specialisation of class ::polar_laurent_coeffs for the ::laplace_3d_HSP_kernel */
template <>
class polar_laurent_coeffs<laplace_3d_HSP_kernel>
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

