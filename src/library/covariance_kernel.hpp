// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2015  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2015  Peter Rucz <rucz@hit.bme.hu>
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
 * \file covariance_kernel.hpp
 * \ingroup library
 * \brief implementation of kernels used to represent covariance function of stochastic processes
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef COVARIANCE_KERNEL_HPP_INCLUDED
#define COVARIANCE_KERNEL_HPP_INCLUDED

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "location_normal.hpp"

namespace NiHu
{

template <class Space>
class covariance_kernel;

/// GENERAL TRAITS
namespace kernel_traits_ns
{
	template <class Space>
	struct space<covariance_kernel<Space> > : Space {};

	template<class Space>
	struct result<covariance_kernel<Space> >
	{
		typedef typename Space::scalar_t type;
	};

	template <class Space>
	struct quadrature_family<covariance_kernel<Space> > : gauss_family_tag {};

	template <class Space>
	struct result_rows<covariance_kernel<Space> > : std::integral_constant<unsigned, 1> {};
	template <class Space>
	struct result_cols<covariance_kernel<Space> > : std::integral_constant<unsigned, 1> {};

	template <class Space>
	struct is_singular<covariance_kernel<Space> > : std::false_type {};

	template <class Space>
	struct test_input<covariance_kernel<Space> > : build<location<Space> > {};

	template <class Space>
	struct trial_input<covariance_kernel<Space> > : build<location<Space> > {};

	template <class Space>
	struct is_symmetric<covariance_kernel<Space> > : std::true_type {};

	/** \todo this is incorrect but works :) */
	template <class Space>
	struct far_field_behaviour<covariance_kernel<Space> > : asymptotic::inverse<1> {};
}


template <class Space>
class covariance_kernel :
	public kernel_base<covariance_kernel<Space> >
{
public:
	typedef kernel_base<covariance_kernel<Space> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t location_t;
	
	covariance_kernel(double sigma, double length) :
		m_sigma(sigma), m_length(length)
	{
	}
	
	static result_t distance(location_t const &x, location_t const &y)
	{
		return (x-y).norm();
	}
	
	static result_t distance_sph(location_t const &x, location_t const &y)
	{
		result_t v = x.normalized().dot(y.normalized());
		v = std::min(v, 1.0);
		v = std::max(v, -1.0);
		return std::acos(v);
	}
	
	result_t operator()(location_t const &x, location_t const &y) const
	{
		auto r = distance_sph(x, y);
		return get_variance() * std::exp(-r / get_correlation_length());
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}
	
	double get_stddev(void) const
	{
		return m_sigma;
	}
	
	double get_variance(void) const
	{
		return m_sigma * m_sigma;
	}
	
	double get_correlation_length(void) const
	{
		return m_length;
	}

private:
	double m_sigma;		/**< \brief Standard deviation */
	double m_length;	/**< \brief Correlation length */
};

} // end of namespace NiHu

#endif // COVARIANCE_KERNEL_HPP_INCLUDED