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
	// the space is inherited
	template <class DK, int Nx, int Ny>
	struct space<normal_derivative_kernel<DK, Nx, Ny> > : space<DK> {};

	// the result is inherited
	template <class DK, int Nx, int Ny>
	struct result<normal_derivative_kernel<DK, Nx, Ny> > : result<DK> {};
	// the number of result rows is inherited
	template <class DK, int Nx, int Ny>
	struct result_rows<normal_derivative_kernel<DK, Nx, Ny> > : result_rows<DK> {};
	// the number of result columns is inherited
	template <class DK, int Nx, int Ny>
	struct result_cols<normal_derivative_kernel<DK, Nx, Ny> > : result_cols<DK> {};

	// the quadrature family is inherited
	template <class DK, int Nx, int Ny>
	struct quadrature_family<normal_derivative_kernel<DK, Nx, Ny> > : quadrature_family<DK> {};

	// singularity is inherited
	template <class DK, int Nx, int Ny>
	struct is_singular<normal_derivative_kernel<DK, Nx, Ny> > : is_singular<DK> {};

	// singular quadrature order is inherited
	template <class DK, int Nx, int Ny>
	struct singular_quadrature_order<normal_derivative_kernel<DK, Nx, Ny> >
		: singular_quadrature_order<DK>  {};

	// singular core is normal derivative of the singular core
	template <class DK, int Nx, int Ny>
	struct singular_core<normal_derivative_kernel<DK, Nx, Ny> > {
		typedef normal_derivative_kernel<
			typename singular_core<DK>::type, Nx, Ny
		> type;
	};

	// symmetric if x and y orders are the same
	template <class DK, int Nx, int Ny>
	struct is_symmetric<normal_derivative_kernel<DK, Nx, Ny> >
		: std::integral_constant<bool, Nx == Ny> {};

	// in the general case, the test input is normal_jacobian
	template <class DK, int Nx, int Ny>
	struct test_input<normal_derivative_kernel<DK, Nx, Ny> >
		: build<location<typename space<DK>::type >, normal_jacobian<typename space<DK>::type > > {};

	// for the Nx = 0 case the test input is location
	template <class DK, int Ny>
	struct test_input<normal_derivative_kernel<DK, 0, Ny> >
		: build<location<typename space<DK>::type> > {};

	// in the general case, the trial input is normal_jacobian
	template <class DK, int Nx, int Ny>
	struct trial_input<normal_derivative_kernel<DK, Nx, Ny> >
		: build<location<typename space<DK>::type >, normal_jacobian<typename space<DK>::type > > {};

	// for the Ny = 0 case the trial input is location
	template <class DK, int Nx>
	struct trial_input<normal_derivative_kernel<DK, Nx, 0> >
		: build<location<typename space<DK>::type> > {};

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
	, public DistanceKernel
{
public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 0, 0> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk = DistanceKernel())
		: DistanceKernel(dk)
	{
	}
	
	result_t operator()(typename base_t::x_t const &x, typename base_t::x_t const &y) const
	{
		x_t rvec = y - x;
		result_t g;
		DistanceKernel::template eval<0>(rvec.norm(), &g);
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
	, public DistanceKernel
{
public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 0, 1> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk = DistanceKernel())
		: DistanceKernel(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &n) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdny = rvec.dot(n) / r;
		result_t f;
		DistanceKernel::template eval<1>(r, &f);
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
	, public DistanceKernel
{
public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 1, 0> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk = DistanceKernel())
		: DistanceKernel(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &n) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdnx = -rvec.dot(n) / r;
		result_t f;
		DistanceKernel::template eval<1>(r, &f);
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
	, public DistanceKernel
{
public:
	typedef kernel_base<normal_derivative_kernel<DistanceKernel, 1, 1> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::scalar_t scalar_t;
	typedef typename base_t::x_t x_t;
	
	normal_derivative_kernel(DistanceKernel const &dk = DistanceKernel())
		: DistanceKernel(dk)
	{
	}
	
	result_t operator()(x_t const &x, x_t const &y, x_t const &nx, x_t const &ny) const
	{
		x_t rvec = y - x;
		scalar_t r = rvec.norm();
		scalar_t rdny = rvec.dot(ny) / r;
		scalar_t rdnx = -rvec.dot(nx) / r;
		result_t f[2];
		DistanceKernel::template eval<2>(r, f);
		return f[0] * rdnx * rdny - f[1] * nx.dot(ny);
	}

	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x(), x.get_unit_normal(), y.get_unit_normal());
	}
};

}


#endif
