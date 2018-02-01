#ifndef NORMAL_DERIVATIVE_KERNEL_HPP_INCLUDED
#define NORMAL_DERIVATIVE_KERNEL_HPP_INCLUDED

#include "../util/brick.hpp"
#include "../core/kernel.hpp"
#include "../library/location_normal.hpp"


namespace NiHu
{
	
/** \brief kernel of the Laplace equation
 * \tparam Space the coordinate space the kernel is defined over
 * \tparam Layer the potential layer tag that can be potential::SLP, potential::DLP, potential::DLPt and potential::HSP
 */
template <class DistanceKernel, int Nx, int Ny>
class normal_derivative_kernel;

/// GENERAL TRAITS
namespace kernel_traits_ns
{
	template <class DK, int Nx, int Ny>
	struct space<normal_derivative_kernel<DK, Nx, Ny> > : space<DK> {};

	template <class DK, int Nx, int Ny>
	struct result<normal_derivative_kernel<DK, Nx, Ny> > : result<DK> {};

	template <class DK, int Nx, int Ny>
	struct quadrature_family<normal_derivative_kernel<DK, Nx, Ny> > : quadrature_family<DK> {};

	template <class DK, int Nx, int Ny>
	struct result_rows<normal_derivative_kernel<DK, Nx, Ny> > : result_rows<DK> {};

	template <class DK, int Nx, int Ny>
	struct result_cols<normal_derivative_kernel<DK, Nx, Ny> > : result_cols<DK> {};

	template <class DK, int Nx, int Ny>
	struct is_singular<normal_derivative_kernel<DK, Nx, Ny> > : is_singular<DK> {};

	template <class DK, int Nx, int Ny>
	struct singular_core<normal_derivative_kernel<DK, Nx, Ny> > {
		typedef normal_derivative_kernel<
			typename singular_core<DK>::type, Nx, Ny
		> type;
	};

	template <class DK, int Nx, int Ny>
	struct is_symmetric<normal_derivative_kernel<DK, Nx, Ny> >
		: std::integral_constant<bool, Nx == Ny> {};

	template <class DK, int Nx, int Ny>
	struct test_input<normal_derivative_kernel<DK, Nx, Ny> > : merge<
		typename test_input<DK>::type,
		typename build<normal_jacobian<typename space<DK>::type> >::type
	> {};

	template <class DK, int Ny>
	struct test_input<normal_derivative_kernel<DK, 0, Ny> >
		: test_input<DK> {};

	template <class DK, int Nx, int Ny>
	struct trial_input<normal_derivative_kernel<DK, Nx, Ny> > : merge<
		typename trial_input<DK>::type,
		typename build<normal_jacobian<typename space<DK>::type> >::type
	> {};

	template <class DK, int Nx>
	struct trial_input<normal_derivative_kernel<DK, Nx, 0> >
		: trial_input<DK> {};

	template <class DK, int Nx, int Ny>
	struct singular_quadrature_order<normal_derivative_kernel<DK, Nx, Ny> >
		: singular_quadrature_order<DK>  {};

	template <class DK>
	struct far_field_behaviour<normal_derivative_kernel<DK, 0, 0> >
		: far_field_behaviour<DK> {};

	template <class DK>
	struct singularity_type<normal_derivative_kernel<DK, 0, 0> >
		: singularity_type<DK> {};
}


template <class DistanceKernel>
class normal_derivative_kernel<DistanceKernel, 0, 0>
	: public kernel_base<normal_derivative_kernel<DistanceKernel, 0, 0> >
{
private:
	DistanceKernel dk;

public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 0, 0> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk) : dk(dk)
	{
	}
	
	result_t operator()(typename base_t::x_t const &x, typename base_t::x_t const &y) const
	{
		x_t rvec = y - x;
		result_t g;
		this->dk.template eval<0>(rvec.norm(), &g);
		return g;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}
};


template <class DistanceKernel>
class normal_derivative_kernel<DistanceKernel, 0, 1>
	: public kernel_base<normal_derivative_kernel<DistanceKernel, 0, 1> >
{
private:
	DistanceKernel dk;

public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 0, 1> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk) : dk(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &n) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdny = rvec.dot(n) / r;
		result_t f;
		this->dk.template eval<1>(r, &f);
		return f * rdny;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x(), y.get_unit_normal());
	}
};

template <class DistanceKernel>
class normal_derivative_kernel<DistanceKernel, 1, 0>
	: public kernel_base<normal_derivative_kernel<DistanceKernel, 1, 0> >
{
private:
	DistanceKernel dk;

public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 1, 0> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk) : dk(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &n) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdnx = -rvec.dot(n) / r;
		result_t f;
		this->dk.template eval<1>(r, &f);
		return f * rdnx;
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x(), x.get_unit_normal());
	}
};


template <class DistanceKernel>
class normal_derivative_kernel<DistanceKernel, 1, 1>
	: public kernel_base<normal_derivative_kernel<DistanceKernel, 1, 1> >
{
private:
	DistanceKernel dk;

public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 1, 1> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk) : dk(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &nx, x_t const &ny) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdny = rvec.dot(ny) / r;
		scalar_t rdnx = -rvec.dot(nx) / r;
		result_t f[2];
		this->dk.template eval<2>(r, f);
		return f[0] * rdnx * rdny - f[1] * nx.dot(ny);
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x(), x.get_unit_normal());
	}
};

}

#include "laplace_base_kernel.hpp"

namespace NiHu
{

/// Laplace Helper Behaviour
namespace kernel_traits_ns
{
	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 0, 1>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 0, 1>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 0, 1>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 0, 1>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 1, 0>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 1, 0>
	> : asymptotic::inverse<1> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 1, 0>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 1, 0>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 1, 1>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_2d<Scalar> >, 1, 1>
	> : asymptotic::inverse<2> {};

	template <class Scalar>
	struct far_field_behaviour<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 1, 1>
	> : asymptotic::inverse<3> {};

	template <class Scalar>
	struct singularity_type<
		normal_derivative_kernel<laplace_helper<space_3d<Scalar> >, 1, 1>
	> : asymptotic::inverse<3> {};
}

}

#endif
