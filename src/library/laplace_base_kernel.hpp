#ifndef LAPLACE_BASE_KERNEL_HPP_INCLUDED
#define LAPLACE_BASE_KERNEL_HPP_INCLUDED

#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "location_normal.hpp"

namespace NiHu
{

template <class Space>
class laplace_helper;

/// GENERAL TRAITS
namespace kernel_traits_ns
{
	template <class Space>
	struct space<laplace_helper<Space> > : Space {};

	template <class Space>
	struct result<laplace_helper<Space> >
	{
		typedef typename Space::scalar_t type;
	};

	template <class Space>
	struct quadrature_family<laplace_helper<Space> > : gauss_family_tag {};

	template <class Space>
	struct result_rows<laplace_helper<Space> > : std::integral_constant<unsigned, 1> {};
	
	template <class Space>
	struct result_cols<laplace_helper<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct is_singular<laplace_helper<Space> > : std::true_type {};

	template <class Space>
	struct singular_core<laplace_helper<Space> > {
		typedef  laplace_helper<Space>  type;
	};

	template <class Space>
	struct test_input<laplace_helper<Space> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<laplace_helper<Space> > : build<location<Space> > {};

	template <class Space>
	struct is_symmetric<laplace_helper<Space> > : std::true_type {};

	template <class Space>
	struct singular_quadrature_order<laplace_helper<Space> >
		: std::integral_constant<unsigned, 7> {};

	template <class Scalar>
	struct far_field_behaviour<laplace_helper<space_2d<Scalar> > >
		: asymptotic::log<1> {};

	template <class Scalar>
	struct singularity_type<laplace_helper<space_2d<Scalar> > >
		: asymptotic::log<1> {};
}

template <class scalar>
class laplace_helper<space_2d<scalar> >
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
		f[1] = 2.0 / (r*r*r) / (2. * M_PI);
		f[0] = -4.0 * f[1];
	}
	
public:
	template <unsigned order>
	void eval(scalar r, scalar *f) const
	{
		this->eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};
	
	
template <class scalar>
class laplace_helper<space_3d<scalar> >
{
	void eval_impl(std::integral_constant<unsigned, 0>, scalar r, scalar *f) const
	{
		*f = 1.0 / r / (4.0 * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 1>, scalar r, scalar *f) const
	{
		*f = -1.0 / r/r / (4.0 * M_PI);
	}
	
	void eval_impl(std::integral_constant<unsigned, 2>, scalar r, scalar *f) const
	{
		f[1] = -1.0 / r/r/r / (4.0 * M_PI);
		f[0] = -3.0 * f[1];
	}
	
	void eval_impl(std::integral_constant<unsigned, 3>, scalar r, scalar *f) const
	{
		f[1] = 3.0 / r/r/r/r / (4.0 * M_PI);
		f[0] = -5 * f[1];
	}
	
public:
	template <unsigned order>
	void eval(scalar r, scalar *f) const
	{
		this->eval_impl(std::integral_constant<unsigned, order>(), r, f);
	}
};

}

#endif
