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

/** \brief the properties of the covariance kernel */
template <class Space>
struct kernel_traits<covariance_kernel<Space> >
{
	typedef typename build<location<Space> >::type test_input_t;
	typedef typename build<location<Space> >::type trial_input_t;
	typedef double result_t;
	enum { result_rows = 1, result_cols = 1 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	
	/** \todo this is incorrect but */
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	
	static bool const is_singular = false;
};

template <class Space>
class covariance_kernel :
	public kernel_base<covariance_kernel<Space> >
{
public:
	typedef kernel_base<covariance_kernel<Space> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::result_t result_t;
	
	covariance_kernel(double sigma, double length) :
		m_sigma(sigma), m_length(length)
	{
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		return get_variance() * std::exp(-r/get_correlation_length());
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

}

#endif // COVARIANCE_KERNEL_HPP_INCLUDED
