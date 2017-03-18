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
#include "potential_kernel.hpp"
#include "location_normal.hpp"

namespace NiHu
{
/** \brief kernel of the Laplace equation
 * \tparam Space the coordinate space the kernel is defined over
 * \tparam Layer the potential layer tag that can be potential::SLP, potential::DLP, potential::DLPt and potential::HSP
 */
template <class Space, class Layer>
class laplace_kernel;


/// GENERAL TRAITS
namespace kernel_traits_ns
{
	template <class Space, class Layer>
	struct space<laplace_kernel<Space, Layer> > : Space {};

	template <class Space, class Layer>
	struct result<laplace_kernel<Space, Layer> >
	{
		typedef typename Space::scalar_t type;
	};

	template <class Space, class Layer>
	struct quadrature_family<laplace_kernel<Space, Layer> > : gauss_family_tag {};

	template <class Space, class Layer>
	struct result_rows<laplace_kernel<Space, Layer> > : std::integral_constant<unsigned, 1> {};
	template <class Space, class Layer>
	struct result_cols<laplace_kernel<Space, Layer> > : std::integral_constant<unsigned, 1> {};

	template <class Space, class Layer>
	struct is_singular<laplace_kernel<Space, Layer> > : std::true_type {};

	template <class Space, class Layer>
	struct singular_core<laplace_kernel<Space, Layer> > {
		typedef  laplace_kernel<Space, Layer>  type;
	};
}

/// SLP TRAITS
namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::SLP> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::SLP> > : build<location<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::SLP> > : std::true_type {};

	/** \brief the singular quadrature order of the laplace SLP kernel
	 * \todo check if the same value can be used for the 2D and 3D case
	 */
	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::SLP> > : std::integral_constant<unsigned, 7> {};
}

/// 2D SLP
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::SLP> > : asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::SLP> > : asymptotic::log<1> {};
}


/** \brief the 2D SLP kernel of the Laplace equation: \f$ -\ln r / 2\pi \f$
 * \tparam Scalar the scalar type
 */
template <class Scalar>
class laplace_kernel<space_2d<Scalar>, potential::SLP>
	: public kernel_base<laplace_kernel<space_2d<Scalar>, potential::SLP> >
{
public:
	typedef kernel_base<laplace_kernel<space_2d<Scalar>, potential::SLP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		Scalar r = (y.get_x() - x.get_x()).norm();
		return -std::log(r) / (2. * M_PI);
	}
};


/// 3D SLP TRAITS
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::SLP> > : asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::SLP> > : asymptotic::inverse<1> {};
}

/** \brief the 3D SLP kernel of the Laplace equation: \f$ 1 / 4 \pi r \f$
 * \tparam Scalar the scalar type
 */
template <class Scalar>
class laplace_kernel<space_3d<Scalar>, potential::SLP>
	: public kernel_base<laplace_kernel<space_3d<Scalar>, potential::SLP> >
{
public:
	typedef kernel_base<laplace_kernel<space_3d<Scalar>, potential::SLP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		Scalar r = (y.get_x() - x.get_x()).norm();
		return 1. / r / (4. * M_PI);
	}
};

/// GENERAL DLP
namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::DLP> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::DLP> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::DLP> > : std::false_type {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::DLP> > : std::integral_constant<unsigned, 7> {};
}

/// 2D DLP
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::DLP> > : asymptotic::inverse<1> {};

	/** \brief kernel singularity type
	 * \todo check this!
	 */
	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::DLP> > : asymptotic::log<1> {};
}

template <class Scalar>
class laplace_kernel<space_2d<Scalar>, potential::DLP>
	: public kernel_base<laplace_kernel<space_2d<Scalar>, potential::DLP> >
{
public:
	typedef kernel_base<laplace_kernel<space_2d<Scalar>, potential::DLP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		typename base_t::x_t rvec = y.get_x() - x.get_x();
		Scalar r = rvec.norm();
		Scalar rdny = rvec.dot(y.get_unit_normal()) / r;
		return -rdny / r / (2. * M_PI);
	}
};

/// 3D DLP
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::DLP> > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::DLP> > : asymptotic::inverse<1> {};
}


template <class Scalar>
class laplace_kernel<space_3d<Scalar>, potential::DLP>
	: public kernel_base<laplace_kernel<space_3d<Scalar>, potential::DLP> >
{
public:
	typedef kernel_base<laplace_kernel<space_3d<Scalar>, potential::DLP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		typename base_t::x_t rvec = y.get_x() - x.get_x();
		Scalar r = rvec.norm();
		auto g = 1. / r / (4. * M_PI);
		Scalar rdny = rvec.dot(y.get_unit_normal()) / r;
		return -g * rdny / r;
	}
};

/// GENERAL DLPt
namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::DLPt> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::DLPt> > : build<location<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::DLPt> > : std::false_type {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::DLPt> > : std::integral_constant<unsigned, 7> {};
}

/// 2D DLPt

namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::DLPt> > : asymptotic::inverse<1> {};

	/** \brief the singularity type
	 * \todo check this just like the plain DLP kernel
	 */
	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::DLPt> > : asymptotic::log<1> {};
}

template <class Scalar>
class laplace_kernel<space_2d<Scalar>, potential::DLPt>
	: public kernel_base<laplace_kernel<space_2d<Scalar>, potential::DLPt> >
{
public:
	typedef kernel_base<laplace_kernel<space_2d<Scalar>, potential::DLPt> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		typename base_t::x_t rvec = trial_input.get_x() - test_input.get_x();
		Scalar r = rvec.norm();
		Scalar rdnx = -rvec.dot(test_input.get_unit_normal()) / r;
		return -rdnx / r / (2. * M_PI);
	}
};

/// 3D DLPt

namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::DLPt> > : asymptotic::inverse<2> {};
	
	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::DLPt> > : asymptotic::inverse<1> {};
}

template <class Scalar>
class laplace_kernel<space_3d<Scalar>, potential::DLPt>
	: public kernel_base<laplace_kernel<space_3d<Scalar>, potential::DLPt> >
{
public:
	typedef kernel_base<laplace_kernel<space_3d<Scalar>, potential::DLPt> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		typename base_t::x_t rvec = trial_input.get_x() - test_input.get_x();
		Scalar r = rvec.norm();
		auto g = 1. / r / (4. * M_PI);
		Scalar rdnx = -rvec.dot(test_input.get_unit_normal()) / r;
		return -g * rdnx / r;
	}
};

/// GENERAL HSP
namespace kernel_traits_ns
{
	template <class Space>
	struct test_input<laplace_kernel<Space, potential::HSP> > : build<location<Space>, normal_jacobian<Space > > {};

	template <class Space>
	struct trial_input<laplace_kernel<Space, potential::HSP> > : build<location<Space>, normal_jacobian<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_kernel<Space, potential::HSP> > : std::true_type {};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space, potential::HSP> > : std::integral_constant<unsigned, 7> {};
}

/// 2D HSP
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar>, potential::HSP> > : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar>, potential::HSP> > : asymptotic::inverse<2> {};
}

template <class Scalar>
class laplace_kernel<space_2d<Scalar>, potential::HSP>
	: public kernel_base<laplace_kernel<space_2d<Scalar>, potential::HSP> >
{
public:
	typedef kernel_base<laplace_kernel<space_2d<Scalar>, potential::HSP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		typename base_t::x_t rvec = trial_input.get_x() - test_input.get_x();
		Scalar r = rvec.norm();
		Scalar rdny = rvec.dot(trial_input.get_unit_normal()) / r;
		Scalar rdnx = -rvec.dot(test_input.get_unit_normal()) / r;
		return 	1./(2. * M_PI)/r/r * (
		test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
		2. * rdny*rdnx);
	}
};

/// 3D HSP

namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar>, potential::HSP> > : asymptotic::inverse<3> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar>, potential::HSP> > : asymptotic::inverse<3> {};
}

template <class Scalar>
class laplace_kernel<space_3d<Scalar>, potential::HSP>
	: public kernel_base<laplace_kernel<space_3d<Scalar>, potential::HSP> >
{
public:
	typedef kernel_base<laplace_kernel<space_3d<Scalar>, potential::HSP> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;

	result_t operator()(test_input_t const &test_input, trial_input_t const &trial_input) const
	{
		typename base_t::x_t rvec = trial_input.get_x() - test_input.get_x();
		Scalar r = rvec.norm();
		Scalar rdny = rvec.dot(trial_input.get_unit_normal()) / r;
		Scalar rdnx = -rvec.dot(test_input.get_unit_normal()) / r;
		auto g = 1. / r / (4. * M_PI);
		return 	g / r / r *
		(
		test_input.get_unit_normal().dot(trial_input.get_unit_normal()) +
		3. * rdnx * rdny
		);
	}
};


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

} // end of namespace NiHu

#include "guiggiani_1992.hpp"

namespace NiHu
{

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

} // end of namespace NiHu

#endif // LAPLACE_KERNEL_HPP_INCLUDED

