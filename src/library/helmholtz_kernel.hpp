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
#include "../util/math_functions.hpp"

#include "laplace_kernel.hpp"

namespace NiHu
{
template <class Space, class WaveNumber>
class helmholtz_helper;

template <class Scalar, class WaveNumber>
class helmholtz_helper<space_2d<Scalar>, WaveNumber>
{
	typedef typename std::complex<Scalar> result_t;
public:
	static void eval(int n, Scalar r, WaveNumber const &k,  result_t *f)
	{
		auto z = k * r;
		auto H0 = bessel::H<0, 2>(std::complex<Scalar>(z));
		f[0] = result_t(0., -1./4.) * H0;
		if (n == 0) return;
		auto H1 = bessel::H<1, 2>(std::complex<Scalar>(z));
		f[1] = k * result_t(0., 1./4.) * H1;
		if (n == 1) return;
		auto H2 = bessel::H<2, 2>(std::complex<Scalar>(z));
		f[2] = k * k * result_t(0., -1./8.) * (H2 - H0);
		if (n == 2) return;
		auto H3 = bessel::H<3, 2>(std::complex<Scalar>(z));
		f[3] = k * k * k * result_t(0., 1./16.) * (H3 - 3. * H1);
	}
};
	
	
template <class Scalar, class WaveNumber>
class helmholtz_helper<space_3d<Scalar>, WaveNumber>
{
public:
	static void eval(int n, Scalar r, WaveNumber const &k, std::complex<Scalar> *f)
	{
		auto kr = k * r;
		auto ikr = std::complex<Scalar>(0,1) * kr;
		f[0] = std::exp(-ikr) / r / (4. * M_PI);
		if (n == 0) return;
		f[1] = f[0] * (-ikr - 1.) / r;
		if (n == 1) return;
		f[2] = f[0] * (2. + 2. * ikr - kr * kr) / r / r;
		if (n == 2) return;
		f[3] = 0;
	}
};
	

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



/// GENERAL DEFINITIONS VALID FOR ALL HELMHOLTZ KERNELS
namespace kernel_traits_ns
{
	// metafunction returning the kernel's space
	template <class Space, class Layer, class WaveNumber>
	struct space<helmholtz_kernel<Space, Layer, WaveNumber> > : Space {};

	// metafunction returning the kernel's result value type
	template <class Space, class Layer, class WaveNumber>
	struct result<helmholtz_kernel<Space, Layer, WaveNumber> >
	{
		typedef std::complex<typename Space::scalar_t> type;
	};

	// metafunction returning the number of rows in the  kernel's result
	template <class Space, class Layer, class WaveNumber>
	struct result_rows<helmholtz_kernel<Space, Layer, WaveNumber> > : std::integral_constant<unsigned, 1> {};

	// metafunction returning the number of columns in the  kernel's result
	template <class Space, class Layer, class WaveNumber>
	struct result_cols<helmholtz_kernel<Space, Layer, WaveNumber> > : std::integral_constant<unsigned, 1> {};

	// metafunction returning the quadrature family the kernel is integrated with
	template <class Space, class Layer, class WaveNumber>
	struct quadrature_family<helmholtz_kernel<Space, Layer, WaveNumber> > : gauss_family_tag {};

	// metafunction indicating if the kernel is singular
	template <class Space, class Layer, class WaveNumber>
	struct is_singular<helmholtz_kernel<Space, Layer, WaveNumber> > : std::true_type {};
	
	// metafunction returning the singular core of the kernel
	template <class Space, class Layer, class WaveNumber>
	struct singular_core<helmholtz_kernel<Space, Layer, WaveNumber> > {
		typedef laplace_kernel<Space, Layer> type;
	};
} // end of namespace kernel_traits_ns



// GENERAL DEFINITIONS FOR THE SLP KERNELS INDEPENDENT OF THE DIMENSION
namespace kernel_traits_ns
{
	// the test input is location
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	// the trial input is location
	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : build<location<Space> > {};

	// the kernel is symmetric
	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : std::true_type {};

	// the singular integrals are performed with 7-th order quadrature
	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::SLP, WaveNumber> > : std::integral_constant<unsigned, 7> {};
}

// SPECIFIC DEFINITIONS FOR THE 2D SLP KERNEL
namespace kernel_traits_ns
{
	// far field behaviour is log r
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::log<1> {};

	// the singularity type is log r
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::log<1> {};
}

// SPECIFIC DEFINITIONS FOR THE 3D SLP KERNEL
namespace kernel_traits_ns
{
	// far field behaviour is 1/r
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::inverse<1> {};

	// the singularity type is 1/r
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::SLP, WaveNumber> > : asymptotic::inverse<1> {};
}


/** \brief Helmholtz SLP kernel
 * \tparam Space the kernel's space type
 * \tparam WaveNumber the wave number's type
 */
template <class Space, class WaveNumber>
class helmholtz_kernel<Space, potential::SLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<Space, potential::SLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
	typedef kernel_base<helmholtz_kernel<Space, potential::SLP, WaveNumber> > base_t;
	
public:
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t x_t;
	
	helmholtz_kernel(WaveNumber const &wave_number)
		: wave_number_kernel<WaveNumber>(wave_number)
	{
	}
	
	/** \brief compute kernel's value from two locations
	 * \param [in] x test location
	 * \param [in] y trial location
	 * \return G(x,y)
	 */
	result_t operator()(x_t const &x, x_t const &y) const
	{
		auto r = (y - x).norm();
		auto const &k = this->get_wave_number();
		result_t g;
		helmholtz_helper<Space, WaveNumber>::eval(0, r, k, &g);
		return g;
	}
	
	/** \brief compute kernel's value from two inputs
	 * \param [in] x test input
	 * \param [in] y trial input
	 * \return G(x,y)
	 */
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}
};



// GENERAL DEFINITIONS FOR THE DLP KERNELS INDEPENDENT OF THE DIMENSION
namespace kernel_traits_ns
{
	// test input is location
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : build<location<Space> > {};

	// trial input is location and normal
	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	// kernel is not symmetric
	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : std::false_type {};

	// the singular integrals are performed with 7-th order quadrature
	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::DLP, WaveNumber> > : std::integral_constant<unsigned, 7> {};
}


// SPECIFIC DEFINITITIONS FOR THE 2D DLP KERNEL
namespace kernel_traits_ns
{
	// far field behaviour is 1/r
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<1> {};

	// singular behaviour is o(1) constant on a smooth surface
	// singular behaviour is o(1/r) on a corner
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<1> {};
}

// SPECIFIC DEFINITIONS FOR THE 3D DLP KERNEL
namespace kernel_traits_ns
{
	// far field behaviour is 1/r^2
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<2> {};

	// singular behaviour is o(1) constant on smooth surfaces
	// singular behaviour is o(1/r^2) on a corner
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::DLP, WaveNumber> > : asymptotic::inverse<2> {};
}



/** \brief kernel of the Helmholtz equation
 * \tparam Space the space type
 * \tparam WaveNumber the wave number type
 */
template <class Space, class WaveNumber>
class helmholtz_kernel<Space, potential::DLP, WaveNumber>
	: public kernel_base<helmholtz_kernel<Space, potential::DLP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<Space, potential::DLP, WaveNumber> > base_t;
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
		auto const &k = this->get_wave_number();
		
		result_t f[2];
		helmholtz_helper<Space, WaveNumber>::eval(1, r, k, f);
		return f[1] * rdny;
	}
};



/// GENERAL DLPt
namespace kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : build<location<Space> > {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : std::false_type {};

	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > : std::integral_constant<unsigned, 7> {};
}


/// 2D DLPt
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<1> {};

	// singular behaviour is o(1) constant on a smooth surface
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::log<1> {};
}


/// 3D DLPt
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<2> {};
	
	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::DLPt, WaveNumber> > : asymptotic::inverse<1> {};
}


/** \brief kernel of the Helmholtz equation
 * \tparam Space the space type
 * \tparam WaveNumber the wave number type
 */
template <class Space, class WaveNumber>
class helmholtz_kernel<Space, potential::DLPt, WaveNumber>
	: public kernel_base<helmholtz_kernel<Space, potential::DLPt, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<Space, potential::DLPt, WaveNumber> > base_t;
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
		auto r = rvec.norm();
		auto rdnx = -rvec.dot(x.get_unit_normal()) / r;
		auto const &k = this->get_wave_number();
		
		result_t f[2];
		helmholtz_helper<Space, WaveNumber>::eval(1, r, k, f);
		return f[1] * rdnx;
	}
};


/// GENERAL HSP
namespace kernel_traits_ns
{
	template <class Space, class WaveNumber>
	struct test_input<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct trial_input<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : build<location<Space>, normal_jacobian<Space>  > {};

	template <class Space, class WaveNumber>
	struct is_symmetric<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : std::true_type {};
	
	template <class Space, class WaveNumber>
	struct singular_quadrature_order<helmholtz_kernel<Space, potential::HSP, WaveNumber> > : std::integral_constant<unsigned, 9> {};
}


/// 2D HSP
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_2d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<2> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_2d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<2> {};
}


/// 3D HSP
namespace kernel_traits_ns
{
	template <class Scalar, class WaveNumber>
	struct far_field_behaviour<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<3> {};

	template <class Scalar, class WaveNumber>
	struct singularity_type<helmholtz_kernel<space_3d<Scalar>, potential::HSP, WaveNumber> > : asymptotic::inverse<3> {};
}



template <class Space, class WaveNumber>
class helmholtz_kernel<Space, potential::HSP, WaveNumber>
	: public kernel_base<helmholtz_kernel<Space, potential::HSP, WaveNumber> >
	, public wave_number_kernel<WaveNumber>
{
public:
	typedef kernel_base<helmholtz_kernel<Space, potential::HSP, WaveNumber> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename Space::scalar_t scalar_t;
	
	helmholtz_kernel(WaveNumber const &k)
		: wave_number_kernel<WaveNumber>(k)
	{
	}

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		auto const &k = this->get_wave_number();
		auto rvec = trial_input.get_x() - test_input.get_x();
		auto r = rvec.norm();
		auto const &nx = test_input.get_unit_normal();
		auto const &ny = trial_input.get_unit_normal();
		scalar_t rdny = rvec.dot(ny) / r;
		scalar_t rdnx = -rvec.dot(nx) / r;

		result_t f[3];
		helmholtz_helper<Space, WaveNumber>::eval(2, r, k, f);
		result_t g1 = f[1] / r;
		result_t g2 = f[2] - g1;
		return 	g2 * rdny*rdnx - g1 * nx.dot(ny);
	}
};



/// USEFUL SHORTHANDS

/** \brief shorthand for the 2d Helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_2d_SLP_kernel = helmholtz_kernel<space_2d<>, potential::SLP, WaveNumber>;
/** \brief shorthand for the 3d Helmholtz SLP kernel */
template <class WaveNumber>
using helmholtz_3d_SLP_kernel = helmholtz_kernel<space_3d<>, potential::SLP, WaveNumber>;
/** \brief shorthand for the 2d Helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_2d_DLP_kernel = helmholtz_kernel<space_2d<>, potential::DLP, WaveNumber>;
/** \brief shorthand for the 3d Helmholtz DLP kernel */
template <class WaveNumber>
using helmholtz_3d_DLP_kernel = helmholtz_kernel<space_3d<>, potential::DLP, WaveNumber>;
/** \brief shorthand for the 2d Helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_2d_DLPt_kernel = helmholtz_kernel<space_2d<>, potential::DLPt, WaveNumber>;
/** \brief shorthand for the 3d Helmholtz DLPt kernel */
template <class WaveNumber>
using helmholtz_3d_DLPt_kernel = helmholtz_kernel<space_3d<>, potential::DLPt, WaveNumber>;
/** \brief shorthand for the 2d Helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_2d_HSP_kernel = helmholtz_kernel<space_2d<>, potential::HSP, WaveNumber>;
/** \brief shorthand for the 3d Helmholtz HSP kernel */
template <class WaveNumber>
using helmholtz_3d_HSP_kernel = helmholtz_kernel<space_3d<>, potential::HSP, WaveNumber>;

} 	// end of namespace NiHu

#endif // HELMHOLTZ_KERNEL_HPP_INCLUDED

