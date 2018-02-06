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
#include "../core/gaussian_quadrature.hpp"
#include "normal_derivative_kernel.hpp"
#include "distance_dependent_kernel.hpp"

namespace NiHu
{
	
template <class Space>
class laplace_kernel;

namespace distance_dependent_kernel_traits_ns
{
	template <class Space>
	struct space<laplace_kernel<Space> > : Space {};

	template <class Space>
	struct result<laplace_kernel<Space> >
	{
		typedef typename Space::scalar_t type;
	};

	template <class Space>
	struct quadrature_family<laplace_kernel<Space> > : gauss_family_tag {};

	template <class Space>
	struct result_rows<laplace_kernel<Space> > : std::integral_constant<unsigned, 1> {};
	
	template <class Space>
	struct result_cols<laplace_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct is_singular<laplace_kernel<Space> > : std::true_type {};

	template <class Space>
	struct singular_core<laplace_kernel<Space> > {
		typedef  laplace_kernel<Space>  type;
	};

	template <class Space>
	struct singular_quadrature_order<laplace_kernel<Space> >
		: std::integral_constant<unsigned, 7> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_2d<Scalar> > >
		: asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_2d<Scalar> > >
		: asymptotic::log<1> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_kernel<space_3d<Scalar> > >
		: asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<laplace_kernel<space_3d<Scalar> > >
		: asymptotic::inverse<1> {};
} // end of distance_dependent_kernel_traits_ns

template <class scalar>
class laplace_kernel<space_2d<scalar> >
	: public distance_dependent_kernel<laplace_kernel<space_2d<scalar> > >
{
private:
	void eval_impl(std::integral_constant<unsigned, 0>, scalar r, scalar *f) const
	{
		*f = -std::log(r) / (2. * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 1>, scalar r, scalar *f) const
	{
		*f = -1.0 / r / (2. * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 2>, scalar r, scalar *f) const
	{
		f[1] = -1.0 / (r*r) / (2. * M_PI);
		f[0] = -2.0  * f[1];
	}
	
	void eval_impl(std::integral_constant<unsigned, 3>, scalar r, scalar *f) const
	{
		auto g = 1. / (r*r*r) / (2. * M_PI);
		f[1] = 2.0 * g; 
		f[0] = -8.0 * g;
	}
	
public:
	template <unsigned order>
	void eval(scalar r, scalar *f) const
	{
		this->eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};
	
	
template <class scalar>
class laplace_kernel<space_3d<scalar> >
	: public distance_dependent_kernel<laplace_kernel<space_3d<scalar> > >
{
	void eval_impl(std::integral_constant<unsigned, 0>, scalar r, scalar *f) const
	{
		*f = 1. / r / (4. * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 1>, scalar r, scalar *f) const
	{
		*f = -1. / (r*r) / (4. * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 2>, scalar r, scalar *f) const
	{
		auto g = 1./(r*r*r)/(4.*M_PI);
		f[1] = -1. * g;
		f[0] = 3. * g;
	}
	
	void eval_impl(std::integral_constant<unsigned, 3>, scalar r, scalar *f) const
	{
		auto g = 1. / (r*r*r*r) / (4. * M_PI);
		f[1] = 3. * g;
		f[0] = -15. * g;
	}
	
public:
	template <unsigned order>
	void eval(scalar r, scalar *f) const
	{
		this->eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};


/// Laplace Helper Behavior
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 0, 1>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 0, 1>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 1, 0>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 1, 0>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 1, 1>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_2d<Scalar> >, 1, 1>
	> : asymptotic::inverse<2> {};
}

namespace kernel_traits_ns
{
	template <class Scalar, int Nx, int Ny>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_kernel<space_3d<Scalar> >, Nx, Ny>
	> : asymptotic::inverse<1 + Nx + Ny> {};

	template <class Scalar, int Nx, int Ny>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_3d<Scalar> >, Nx, Ny>
	> : asymptotic::inverse<1 + Nx + Ny> {};

	// the normal derivative cancels the -2 singularity on a smooth boundary
	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_3d<Scalar> >, 0, 1>
	> : asymptotic::inverse<1> {};

	// the normal derivative cancels the -2 singularity on a smooth boundary
	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_kernel<space_3d<Scalar> >, 1, 0>
	> : asymptotic::inverse<1> {};
}
	

/** \brief shorthand for the 2d Laplace SLP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_2d<> >, 0, 0> laplace_2d_SLP_kernel;
/** \brief shorthand for the 3d Laplace SLP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_3d<> >, 0, 0> laplace_3d_SLP_kernel;
/** \brief shorthand for the 2d Laplace DLP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_2d<> >, 0, 1> laplace_2d_DLP_kernel;
/** \brief shorthand for the 3d Laplace DLP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_3d<> >, 0, 1> laplace_3d_DLP_kernel;
/** \brief shorthand for the 2d Laplace DLPt kernel */
typedef normal_derivative_kernel<laplace_kernel<space_2d<> >, 1, 0> laplace_2d_DLPt_kernel;
/** \brief shorthand for the 3d Laplace DLPt kernel */
typedef normal_derivative_kernel<laplace_kernel<space_3d<> >, 1, 0> laplace_3d_DLPt_kernel;
/** \brief shorthand for the 2d Laplace HSP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_2d<> >, 1, 1> laplace_2d_HSP_kernel;
/** \brief shorthand for the 3d Laplace HSP kernel */
typedef normal_derivative_kernel<laplace_kernel<space_3d<> >, 1, 1> laplace_3d_HSP_kernel;

} // end of namespace NiHu

#include "guiggiani_1992.hpp"

namespace NiHu
{

/** \brief specialisation of class ::polar_laurent_coeffs for the ::laplace_3d_HSP_kernel */
template <class Scalar>
class polar_laurent_coeffs<
	normal_derivative_kernel<laplace_kernel<space_3d<Scalar> >, 1, 1>
>
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

