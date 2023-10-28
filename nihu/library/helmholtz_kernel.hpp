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
 * \brief Kernels of the Helmholtz equation \f$ \nabla^2 p + k^2 p = 0 \f$
 * \ingroup lib_helmholtz
 */

#ifndef HELMHOLTZ_KERNEL_HPP_INCLUDED
#define HELMHOLTZ_KERNEL_HPP_INCLUDED

#include <boost/math/constants/constants.hpp>

#include <cmath>

#include "../core/global_definitions.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/math_functions.hpp"
#include "nihu/library/normal_derivative_kernel.hpp"
#include "distance_dependent_kernel.hpp"

#include "laplace_kernel.hpp"

namespace NiHu
{
template <class Space, class WaveNumber>
class helmholtz_kernel;

/// GENERAL TRAITS
namespace distance_dependent_kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct space<helmholtz_kernel<Space, WaveNumber> > : Space {};

	template <class Space, class WaveNumber>
	struct result<helmholtz_kernel<Space, WaveNumber> >
	{
		typedef std::complex<typename Space::scalar_t> type;
	};

	template <class Space, class WaveNumber>
	struct quadrature_family<helmholtz_kernel<Space, WaveNumber> >
		: gauss_family_tag {};

	template <class Space, class WaveNumber>
	struct is_singular<helmholtz_kernel<Space, WaveNumber> > : std::true_type {};

	template <class Space, class WaveNumber>
	struct singular_core<helmholtz_kernel<Space, WaveNumber> > {
		typedef laplace_kernel<Space>  type;
	};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, WaveNumber> >
		: std::integral_constant<unsigned, 7> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, WaveNumber> >
		: asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, WaveNumber> >
		: asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, WaveNumber> >
		: asymptotic::inverse<1> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, WaveNumber> >
		: asymptotic::inverse<1> {};
} // distance_dependent_kernel_traits_ns


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


template <class scalar, class WaveNumber>
class helmholtz_kernel<space_2d<scalar>, WaveNumber>
	: public wave_number_kernel<WaveNumber>
	, public distance_dependent_kernel<helmholtz_kernel<space_2d<scalar>, WaveNumber> >
{
public:
	typedef typename distance_dependent_kernel_traits_ns::result<
		helmholtz_kernel<space_2d<scalar>, WaveNumber>
	>::type result_t;
	
private:
	void eval_impl(std::integral_constant<unsigned, 0>, scalar r, result_t *f) const
	{
		auto z = this->get_wave_number() * r;
		auto H0 = bessel::H<0, 2>(std::complex<scalar>(z));
		*f = std::complex<scalar>(0, -.25) * H0;
	}
	
	void eval_impl(std::integral_constant<unsigned, 1>, scalar r, result_t *f) const
	{
		auto const &k = this->get_wave_number();
		auto z = k * r;
		auto H1 = bessel::H<1, 2>(std::complex<scalar>(z));
		*f = std::complex<scalar>(0, .25) * this->get_wave_number() * H1;
	}
	
	void eval_impl(std::integral_constant<unsigned, 2>, scalar r, result_t *f) const
	{
		auto const &k = this->get_wave_number();
		auto z = k * r;
		auto H0 = bessel::H<0, 2>(std::complex<scalar>(z));
		auto H1 = bessel::H<1, 2>(std::complex<scalar>(z));
		auto H2 = bessel::H<2, 2>(std::complex<scalar>(z));
		f[1] = std::complex<scalar>(0, .25) * k*k * (H1/z);
		f[0] = std::complex<scalar>(0, -.25) * k*k * (.5 * H2 - .5 * H0 + H1/z);
	}
	
	void eval_impl(std::integral_constant<unsigned, 3>, scalar r, result_t *f) const
	{
		f[0] = f[1] = result_t(0);
	}
	
public:
	helmholtz_kernel(WaveNumber const &k)
		: wave_number_kernel<WaveNumber>(k)
	{
	}
	
	template <unsigned order>
	void eval(scalar r, std::complex<scalar> *f) const
	{
		eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};
	
	
template <class scalar, class WaveNumber>
class helmholtz_kernel<space_3d<scalar>, WaveNumber>
	: public wave_number_kernel<WaveNumber>
	, public distance_dependent_kernel<helmholtz_kernel<space_3d<scalar>, WaveNumber> >
{
public:
	typedef typename distance_dependent_kernel_traits_ns::result<
		helmholtz_kernel<space_3d<scalar>, WaveNumber>
	>::type result_t;
	
private:
	void eval_impl(std::integral_constant<unsigned, 0>, scalar r, result_t *f) const
	{
		using namespace boost::math::double_constants;
		auto kr = this->get_wave_number() * r;
		auto ikr = result_t(0,1) * kr;
		f[0] = std::exp(-ikr) / r / (4. * pi);
	}
	

	void eval_impl(std::integral_constant<unsigned, 1>, scalar r, result_t *f) const
	{
		using boost::math::double_constants::pi;
		auto kr = this->get_wave_number() * r;
		auto ikr = result_t(0,1) * kr;
		f[0] = std::exp(-ikr) / (r*r) / (4. * pi) * (-ikr - 1.);
	}
	

	void eval_impl(std::integral_constant<unsigned, 2>, scalar r, result_t *f) const
	{
		using boost::math::double_constants::pi;
		auto ikr = result_t(0,1) * (this->get_wave_number() * r);
		auto g = std::exp(-ikr) / (r*r*r) / (4. * pi);
		f[1] = -g * (1. + ikr);
		f[0] = g * (3. + ikr * (3. + ikr));
	}
	

	void eval_impl(std::integral_constant<unsigned, 3>, scalar r, result_t *f) const
	{
		using boost::math::double_constants::pi;
		auto ikr = result_t(0,1) * (this->get_wave_number() * r);
		auto g = std::exp(-ikr)/(r*r*r*r) / (4. * pi);
		f[1] = g * (3. + ikr * (3. + ikr));
		f[0] = -g * (15. + ikr * (15 + ikr * (6 + ikr)));
	}
	
public:
	helmholtz_kernel(WaveNumber const &k)
		: wave_number_kernel<WaveNumber>(k)
	{
	}
	
	template <unsigned order>
	void eval(scalar r, result_t *f) const
	{
		this->eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};

/// Helmholtz normal derivative kernel behaviors
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 0, 0>
	> : asymptotic::log<1> {};

	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 0, 0>
	> : asymptotic::log<1> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 0, 1>
	> : asymptotic::inverse<1> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 0, 1>
	> : asymptotic::inverse<1> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 1, 0>
	> : asymptotic::inverse<1> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 1, 0>
	> : asymptotic::inverse<1> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 1, 1>
	> : asymptotic::inverse<2> {};

	/// \todo check this
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_2d<Scalar>, WaveNumber>, 1, 1>
	> : asymptotic::inverse<2> {};
	
}

/// Helmholtz normal derivative kernel behaviors
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber, int Nx, int Ny>
	struct far_field_behaviour<
		normal_derivative_kernel<helmholtz_kernel<space_3d<Scalar>, WaveNumber>, Nx, Ny>
	> : asymptotic::inverse<1 + Nx + Ny> {};

	template <class Scalar, class WaveNumber, int Nx, int Ny>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_3d<Scalar>, WaveNumber>, Nx, Ny>
	> : asymptotic::inverse<1 + Nx + Ny> {};

	// the normal derivative on the smooth boundary cancels the 1/r^2
	// singularity of the kernel's derivative
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_3d<Scalar>, WaveNumber>, 0, 1>
	> : asymptotic::inverse<1> {};

	// the normal derivative on the smooth boundary cancels the 1/r^2
	// singularity of the kernel's derivative
	template <class Scalar, class WaveNumber>
	struct singularity_type<
		normal_derivative_kernel<helmholtz_kernel<space_3d<Scalar>, WaveNumber>, 1, 0>
	> : asymptotic::inverse<1> {};
} // end of namespace kernel_traits_ns



/** \brief shorthand for the 2d Helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_2d_SLP_kernel = normal_derivative_kernel<helmholtz_kernel<space_2d<>, WaveNumber>, 0, 0>;
/** \brief shorthand for the 3d Helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_3d_SLP_kernel = normal_derivative_kernel<helmholtz_kernel<space_3d<>, WaveNumber>, 0, 0>;
/** \brief shorthand for the 2d Helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_2d_DLP_kernel = normal_derivative_kernel<helmholtz_kernel<space_2d<>, WaveNumber>, 0, 1>;
/** \brief shorthand for the 3d Helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_3d_DLP_kernel = normal_derivative_kernel<helmholtz_kernel<space_3d<>, WaveNumber>, 0, 1>;
/** \brief shorthand for the 2d Helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_2d_DLPt_kernel = normal_derivative_kernel<helmholtz_kernel<space_2d<>, WaveNumber>, 1, 0>;
/** \brief shorthand for the 3d Helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_3d_DLPt_kernel = normal_derivative_kernel<helmholtz_kernel<space_3d<>, WaveNumber>, 1, 0>;
/** \brief shorthand for the 2d Helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_2d_HSP_kernel = normal_derivative_kernel<helmholtz_kernel<space_2d<>, WaveNumber>, 1, 1>;
/** \brief shorthand for the 3d Helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_3d_HSP_kernel = normal_derivative_kernel<helmholtz_kernel<space_3d<>, WaveNumber>, 1, 1>;
/** \brief shorthand for the 3d Helmholtz xx kernel */
template <class WaveNumber>
using helmholtz_3d_xx_kernel = normal_derivative_kernel<helmholtz_kernel<space_3d<>, WaveNumber>, 2, 0>;

} 	// end of namespace NiHu

#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED

