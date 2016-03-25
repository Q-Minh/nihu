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
 * \file helmholtz_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels of the Helmholtz equation \f$ \nabla^2 p + k^2 p = 0 \f$
 */

#ifndef HELMHOLTZ_KERNEL_HPP_INCLUDED
#define HELMHOLTZ_KERNEL_HPP_INCLUDED

#include <cmath>
#include <complex>

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"

#include "location_normal.hpp"

#include "basic_bricks.hpp"
#include "../util/collection.hpp"
#include "../util/math_functions.hpp"

#include "laplace_kernel.hpp"

/**
 * \brief kernel data that stores the wave number
 * \tparam wave_number_type the wave number type
 */
template <class wave_number_type>
class wave_number_data
{
public:
	/** \brief constructor setting the wave number
	 * \param [in] wn the wave number to set
	 */
	wave_number_data(wave_number_type const &wn = wave_number_type()) :
		m_wave_number(wn)
	{
	}

	/** \brief return wave number
	 * \return wave number
	 */
	wave_number_type const &get_wave_number(void) const
	{
		return m_wave_number;
	}

private:
	wave_number_type m_wave_number;
};

/** \brief a brick representing the expression \f$ i k r \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct ikr_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_ikr(std::complex<scalar>(0.,1.) * kernel_data.get_wave_number() * wall::get_distance())
		{
		}

		/** \brief return ikr
		 * \return ikr
		 */
		result_t const &get_ikr(void) const
		{
			return m_ikr;
		}

	private:
		result_t m_ikr;
	};
};



/** \brief a brick representing the expression \f$ k r \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct kr_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_kr(kernel_data.get_wave_number() * wall::get_distance())
		{
		}

		/** \brief return kr
		 * \return kr
		 */
		result_t const &get_kr(void) const
		{
			return m_kr;
		}

	private:
		result_t m_kr;
	};
};


/** \brief a brick representing the expression \f$ k r H_0^(2)(kr) \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct krH0_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_krH0(wall::get_kr() * bessel::H<0,2>(wall::get_kr()))
		{
		}

		/** \brief return kr
		 * \return kr
		 */
		result_t const &get_krH0(void) const
		{
			return m_krH0;
		}

	private:
		result_t m_krH0;
	};
};


/** \brief a brick representing the expression \f$ H_1^(2)(kr) \f$
 * \tparam scalar scalar type of the coordinate space the distance is defined over
 */
template <class scalar>
struct H1_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_H1(bessel::H<1,2>(wall::get_kr()))
		{
		}

		/** \brief return kr
		 * \return kr
		 */
		result_t const &get_H1(void) const
		{
			return m_H1;
		}

	private:
		result_t m_H1;
	};
};


/** \brief a brick representing a Helmholtz kernel
 * \tparam Space the coordinate space the distance is defined over
 * \tparam Layer the potential type tag
 */
template <class Space, class Layer>
struct helmholtz_brick;

/** \brief combination of several bricks into a Helmholtz wall
 * \tparam Space the coordinate space the Helmholtz kernel is defined over
 * \tparam Layer the potential type tag
 */
template <class Space, class Layer>
struct helmholtz_wall;

/** \brief kernel of the Helmholtz equation
 * \tparam Space the coordinate space
 * \tparam Layer the potential type tag
 * \tparam WaveNumber the wave number type
 */
template <class Space, class Layer, class WaveNumber>
class helmholtz_kernel :
	public kernel_base<helmholtz_kernel<Space, Layer, WaveNumber> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_kernel(WaveNumber const &wave_number) :
		kernel_base<helmholtz_kernel<Space, Layer, WaveNumber> >(wave_number_data<WaveNumber>(wave_number))
	{
	}
};

/** \brief specialisation of ::helmholtz_brick for the 2D SLP case */
template <class Scalar>
struct helmholtz_brick<space_2d<Scalar>, potential::SLP>
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_g(result_t(0., -.25) * bessel::H<0, 2>(wall::get_kr()))
		{
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_helmholtz_g(void) const
		{
			return m_helmholtz_g;
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_g;
		}

	private:
		result_t m_helmholtz_g;
	};
};


/** \brief specialisation of ::helmholtz_brick for the 3D SLP case */
template <class Scalar>
struct helmholtz_brick<space_3d<Scalar>, potential::SLP>
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_g(std::exp(-wall::get_ikr()) / wall::get_distance() / (4. * M_PI))
		{
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_helmholtz_g(void) const
		{
			return m_helmholtz_g;
		}

		/** \brief return Helmholtz g kernel
		 * \return Helmholtz g kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_g;
		}

	private:
		result_t m_helmholtz_g;
	};
};


/** \brief specialisation of ::helmholtz_wall for the 2D SLP case */
template <class Scalar>
struct helmholtz_wall<space_2d<Scalar>, potential::SLP> : build<
	distance_vector_brick<space_2d<Scalar> >,
	distance_brick<Scalar>,
	kr_brick<Scalar>,
	helmholtz_brick<space_2d<Scalar>, potential::SLP>
> {};

/** \brief specialisation of ::helmholtz_wall for the 3D SLP case */
template <class Scalar>
struct helmholtz_wall<space_3d<Scalar>, potential::SLP> : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	ikr_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::SLP>
> {};

namespace kernel_traits_ns
{
	template <class Space, class Layer, class WaveNumber>
	struct space<helmholtz_kernel<Space, Layer, WaveNumber> > : Space {};

	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	template <class Space, class Layer, class WaveNumber>
	struct data<helmholtz_kernel<Space, Layer, WaveNumber> > {
		typedef collect<wave_number_data<WaveNumber> > type;
	};

	template <class Space, class Layer, class WaveNumber>
	struct output<helmholtz_kernel<Space, Layer, WaveNumber> > : helmholtz_wall<Space, Layer> {};

	template <class Space, class Layer, class WaveNumber>
	struct result_rows<helmholtz_kernel<Space, Layer, WaveNumber> > : std::integral_constant<unsigned, 1> {};
	template <class Space, class Layer, class WaveNumber>
	struct result_cols<helmholtz_kernel<Space, Layer, WaveNumber> > : std::integral_constant<unsigned, 1> {};

	template <class Space, class Layer, class WaveNumber>
	struct quadrature_family<helmholtz_kernel<Space, Layer, WaveNumber> > : gauss_family_tag {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : std::true_type {};

	template <class Space, class Layer, class WaveNumber>
	struct is_singular<helmholtz_kernel<Space, Layer, WaveNumber> > : std::true_type {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : std::integral_constant<unsigned, 7> {};

	template <class Space, class Layer, class WaveNumber>
	struct singular_core<helmholtz_kernel<Space, Layer, WaveNumber> > {
		typedef laplace_kernel<Space, Layer> type;
	};
}


/** \brief specialisation of ::helmholtz_brick for the 2d DLP case */
template <class Scalar>
struct helmholtz_brick<space_2d<Scalar>, potential::DLP>
{
	/** \brief the brick template
	* \tparam the wall the brick is placed on
	*/
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		* \tparam test_input_t the test input type
		* \tparam trial_input_t the trial input type
		* \tparam kernel_data_t the kernel data type
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
			m_helmholtz_h(result_t(0., .25) * kernel_data.get_wave_number() *
			wall::get_H1() * wall::get_rdny())
		{
		}

		/** \brief return Helmholtz h kernel
		* \return Helmholtz h kernel
		*/
		result_t const &get_helmholtz_h(void) const
		{
			return m_helmholtz_h;
		}

		/** \brief return Helmholtz h kernel
		* \return Helmholtz h kernel
		*/
		result_t const &get_result(void) const
		{
			return m_helmholtz_h;
		}

	private:
		result_t m_helmholtz_h;
	};
};


/** \brief specialisation of ::helmholtz_brick for the 3d DLP case */
template <class Scalar>
struct helmholtz_brick<space_3d<Scalar>, potential::DLP>
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_h(-(1.+wall::get_ikr()) * wall::get_helmholtz_g() / wall::get_distance() * wall::get_rdny())
		{
		}

		/** \brief return Helmholtz DLP kernel
		 * \return Helmholtz DLP kernel
		 */
		result_t const &get_helmholtz_h(void) const
		{
			return m_helmholtz_h;
		}

		/** \brief return Helmholtz DLP kernel
		 * \return Helmholtz DLP kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_h;
		}

	private:
		result_t m_helmholtz_h;
	};
};


/** \brief specialisation of ::helmholtz_wall for the 2D DLP case */
template <class Scalar>
struct helmholtz_wall<space_2d<Scalar>, potential::DLP> : build<
	distance_vector_brick<space_2d<Scalar> >,
	distance_brick<Scalar>,
	kr_brick<Scalar>,
	H1_brick<Scalar>,
	rdny_brick<Scalar>,
	helmholtz_brick<space_2d<Scalar>, potential::DLP>
> {};


/** \brief specialisation of ::helmholtz_wall for the 3D DLP case */
template <class Scalar>
struct helmholtz_wall<space_3d<Scalar>, potential::DLP> : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	ikr_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::SLP>,
	rdny_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::DLP>
> {};


namespace kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : build<location<Space> > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : std::false_type {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<2> {};

	/** \todo check this */
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : std::integral_constant<unsigned, 7> {};
}

/** \brief specialisation og ::helmholtz_brick for the 3D DLPt case */
template <class Scalar>
struct helmholtz_brick<space_3d<Scalar>, potential::DLPt>
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_ht(-(1.+wall::get_ikr()) * wall::get_helmholtz_g() / wall::get_distance() * wall::get_rdnx())
		{
		}

		/** \brief return Helmholtz Ht kernel
		 * \return Helmholtz Ht kernel
		 */
		result_t const &get_helmholtz_ht(void) const
		{
			return m_helmholtz_ht;
		}

		/** \brief return Helmholtz h kernel
		 * \return Helmholtz h kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_ht;
		}

	private:
		result_t m_helmholtz_ht;
	};
};

/** \brief specialisation of ::helmholtz_wall for the 3D DLPt case */
template <class Scalar>
struct helmholtz_wall<space_3d<Scalar>, potential::DLPt> : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	ikr_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::SLP>,
	rdnx_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::DLPt>
> {};

namespace kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : build<location<Space> > {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : std::false_type {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<2> {};

	/** \todo check this */
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<1> {};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : std::integral_constant<unsigned, 7> {};
}

/** \brief specialisation of ::helmholtz_brick for the 3D HSP case */
template <class Scalar>
struct helmholtz_brick<space_3d<Scalar>, potential::HSP>
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<Scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_hyper(
				wall::get_helmholtz_g()/wall::get_distance()/wall::get_distance() * (
					(1. + wall::get_ikr())*test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
					(3. + 3.*wall::get_ikr() + wall::get_ikr()*wall::get_ikr())*wall::get_rdny()*wall::get_rdnx()
				)
			)
		{
		}

		/** \brief return Helmholtz hypersingular kernel
		 * \return Helmholtz hypersingular kernel
		 */
		result_t const &get_helmholtz_hyper(void) const
		{
			return m_helmholtz_hyper;
		}

		/** \brief return Helmholtz hypersingular kernel
		 * \return Helmholtz hypersingular kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_hyper;
		}

	private:
		result_t m_helmholtz_hyper;
	};
};

/** \brief specialisation of ::helmholtz_wall for the 3D HSP case */ 
template <class Scalar>
struct helmholtz_wall<space_3d<Scalar>, potential::HSP > : build<
	distance_vector_brick<space_3d<Scalar> >,
	distance_brick<Scalar>,
	ikr_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::SLP>,
	rdnx_brick<Scalar>,
	rdny_brick<Scalar>,
	helmholtz_brick<space_3d<Scalar>, potential::HSP>
> {};


namespace kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : std::true_type {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<2> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<3> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::log<2> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<3> {};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : std::integral_constant<unsigned, 9> {};
}


/** \brief a brick representing a 2D Helmholtz double derivative kernel
 * \tparam scalar the scalar of the coordinate space the distance is defined over
 */
template <class scalar, unsigned i, unsigned j>
struct helmholtz_2d_double_brick
{
	/** \brief the brick template
	 * \tparam the wall the brick is placed on
	 */
	template <class wall>
	class brick : public wall
	{
	public:
		/** \brief the result type */
		typedef std::complex<scalar> result_t;

		/** \brief templated constructor
		 * \tparam test_input_t the test input type
		 * \tparam trial_input_t the trial input type
		 * \tparam kernel_data_t the kernel data type
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
			m_helmholtz_double(
				result_t(0., .25) * kernel_data.get_wave_number() / wall::get_distance() * (
					wall::get_distance_vector()(i) / wall::get_distance() *
					wall::get_distance_vector()(j) / wall::get_distance() *
					( wall::get_krH0() - 2.*wall::get_H1() ) +
					(i == j ? wall::get_H1() : 0.)
				)
			)
		{
		}

		/** \brief return Helmholtz 2d double kernel
		 * \return Helmholtz double kernel
		 */
		result_t const &get_helmholtz_double(void) const
		{
			return m_helmholtz_double;
		}

		/** \brief return Helmholtz 2d double kernel
		 * \return Helmholtz 2d double kernel
		 */
		result_t const &get_result(void) const
		{
			return m_helmholtz_double;
		}

	private:
		result_t m_helmholtz_double;
	};
};


template <class scalar, int i, int j>
struct helmholtz_2d_double_wall : build<
	distance_vector_brick<space_2d<scalar> >,
	distance_brick<scalar>,
	kr_brick<scalar>,
	krH0_brick<scalar>,
	H1_brick<scalar>,
	helmholtz_2d_double_brick<scalar, i, j>
> {};


// forward declaration
template <class wave_number_t, int i, int j>
class helmholtz_2d_double_kernel;

/** \brief traits of the Helmholtz Hyper kernel */
template <class wave_number_t, int i, int j>
struct kernel_traits<helmholtz_2d_double_kernel<wave_number_t, i, j> >
{
	/** \brief kernel test input type */
	typedef location_input_2d test_input_t;
	/** \brief kernel trial input type */
	typedef location_input_2d trial_input_t;
	/** \brief kernel data type */
	typedef collect<wave_number_data<wave_number_t> > data_t;
	/** \brief the kernel output type */
	typedef typename helmholtz_2d_double_wall<space_2d<>::scalar_t, i, j>::type output_t;
	/** \brief the kernel result's dimension */
	enum { result_rows = 1, result_cols = 1 };
	/** \brief the quadrature family the kernel is integrated with */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief indicates if K(x,y) = K(y,x) */
	static bool const is_symmetric = true;
	/** \brief indicates whether kernel is singular */
	static bool const is_singular = true;

	/** \brief the far field asymptotic behaviour of the kernel */
	typedef asymptotic::inverse<2> far_field_behaviour_t;
};


/** \brief singular traits of the Helmholtz 2d double kernel */
template <class wave_number_t, int i, int j>
struct singular_kernel_traits<helmholtz_2d_double_kernel<wave_number_t, i, j> >
{
	/** \brief singularity type */
	typedef asymptotic::inverse<2> singularity_type_t;
	/** \brief the singularity type when used with guiggiani's method */
//	typedef laplace_2d_double_kernel<i, j> singular_core_t;
	/** \brief quadrature order  */
	static unsigned const singular_quadrature_order = 9;
};


/** \brief Double derivative kernel of the Helmholtz equation in 2D */
template <class wave_number_t, int i, int j>
class helmholtz_2d_double_kernel :
	public kernel_base<helmholtz_2d_double_kernel<wave_number_t, i, j> >
{
public:
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_2d_double_kernel(wave_number_t const &wave_number) :
		kernel_base<helmholtz_2d_double_kernel<wave_number_t, i, j> >(wave_number_data<wave_number_t>(wave_number))
	{
	}
};


/** \brief shorthand for the 2d helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_2d_SLP_kernel = helmholtz_kernel<space_2d<>, potential::SLP, WaveNumber>;
/** \brief shorthand for the 3d helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_3d_SLP_kernel = helmholtz_kernel<space_3d<>, potential::SLP, WaveNumber>;
/** \brief shorthand for the 2d helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_2d_DLP_kernel = helmholtz_kernel<space_2d<>, potential::DLP, WaveNumber>;
/** \brief shorthand for the 3d helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_3d_DLP_kernel = helmholtz_kernel<space_3d<>, potential::DLP, WaveNumber>;
/** \brief shorthand for the 2d helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_2d_DLPt_kernel = helmholtz_kernel<space_2d<>, potential::DLPt, WaveNumber>;
/** \brief shorthand for the 3d helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_3d_DLPt_kernel = helmholtz_kernel<space_3d<>, potential::DLPt, WaveNumber>;
/** \brief shorthand for the 2d helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_2d_HSP_kernel = helmholtz_kernel<space_2d<>, potential::HSP, WaveNumber>;
/** \brief shorthand for the 3d helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_3d_HSP_kernel = helmholtz_kernel<space_3d<>, potential::HSP, WaveNumber>;


#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED

