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

namespace NiHu
{

template <class wave_number_type>
class wave_number_kernel
{
public:
	/** \brief constructor setting the wave number
	 * \param [in] wn the wave number to set
	 */
	wave_number_kernel(wave_number_type const &wn)
		: m_wave_number(wn)
	{
	}

	wave_number_type const &get_wave_number(void) const
	{
		return m_wave_number;
	}

private:
	wave_number_type const m_wave_number;
};


/** \brief kernel of the Helmholtz equation
 * \tparam Space the coordinate space
 * \tparam Layer the potential type tag
 * \tparam WaveNumber the wave number type
 */
template <class Space, class Layer, class WaveNumber>
class helmholtz_kernel;


template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto r = (y.get_x() - x.get_x()).norm();
		auto kr = this->get_wave_number() * r;
		return result_t(0., -.25) * bessel::H<0, 2>(result_t(kr));
	}
};


template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto r = (y.get_x() - x.get_x()).norm();
		auto ikr = std::complex<Scalar>(0,1) * (this->get_wave_number() * r);
		return std::exp(-ikr) / r / (4. * M_PI);
	}
};


template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = (rvec).norm();
		auto rdny = rvec.dot(y.get_unit_normal()) / r;
		auto ikr = std::complex<Scalar>(0,1) * (this->get_wave_number() * r);
		auto g = std::exp(-ikr) / r / (4. * M_PI);
		return -(1.+ikr) * g / r * rdny;
	}
};


template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = (rvec).norm();
		auto rdnx = -rvec.dot(x.get_unit_normal()) / r;
		auto ikr = std::complex<Scalar>(0,1) * (this->get_wave_number() * r);
		auto g = std::exp(-ikr) / r / (4. * M_PI);
		return -(1.+ikr) * g / r * rdnx;
	}
};




namespace kernel_traits_ns
{
	template <class Space, class Layer, class WaveNumber>
	struct space<helmholtz_kernel<Space, Layer, WaveNumber> > : Space {};

	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	template <class Space, class Layer, class WaveNumber>
	struct result<helmholtz_kernel<Space, Layer, WaveNumber> >
	{
		typedef std::complex<typename Space::scalar_t> type;
	};

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


/** \brief kernel of the Helmholtz equation
 * \tparam Space the coordinate space
 * \tparam Layer the potential type tag
 * \tparam WaveNumber the wave number type
 */
template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto rdny = rvec.dot(y.get_unit_normal()) / r;
		auto kr = this->get_wave_number() * r;
		auto H1 = bessel::H<1,2>(result_t(kr));
		return result_t(0., .25) * this->get_wave_number() * H1 * rdny;
	}
};

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


template <class Scalar, class WaveNumber>
class helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber>
	: public kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	/** \brief constructor
	 * \param [in] wave_number the wave number
	 */
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		auto ikr = std::complex<Scalar>(0,1) * this->get_wave_number() * r;
		auto rdny = rvec.dot(y.get_unit_normal()) / r;
		auto rdnx = -rvec.dot(x.get_unit_normal()) / r;
		auto g = std::exp(-ikr) / r / (4. * M_PI);
		
		return g / r / r * (
					(1. + ikr)*x.get_unit_normal().dot(y.get_unit_normal()) +
					(3. + 3.*ikr + ikr*ikr)*rdny*rdnx
				);
	}
};

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

} 	// end of namespace NiHu

#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED

