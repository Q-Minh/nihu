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
 * \brief implementation of kernels used to represente covariance function of stochastic processes
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef COVARIANCE_KERNEL_HPP_INCLUDED
#define COVARIANCE_KERNEL_HPP_INCLUDED

#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "../util/collection.hpp"
#include "location_normal.hpp"
#include "basic_bricks.hpp"

namespace NiHu
{

class covariance_data
{
public:
	covariance_data(double sigma, double d)
		: m_sigma(sigma), m_length(d)
	{
	}
	
	double get_correlation_length(void) const
	{
		return m_length;
	}
	
	double get_stddev(void) const
	{
		return m_sigma;
	}
	
	double get_variance(void) const
	{
		return m_sigma * m_sigma;
	}

private:
	double m_sigma;		/**< \brief Standard deviation */
	double m_length;	/**< \brief Correlation length */
};

template <class Space>
struct CKernel
{
	typedef Eigen::Matrix<double, 1, 1> return_type;
	
	typedef typename build<location<Space> >::type location_input;
	
	return_type operator()(
		location_input const &x,
		location_input const &y,
		covariance_data const &data)
	{
		auto rvec = y.get_x() - x.get_x();
		auto r = rvec.norm();
		return_type ret;
		ret(0,0) = data.get_variance() * std::exp(-r/data.get_correlation_length());
		return ret;
	}
};

template <class Space>
class covariance_kernel;

/** \brief the properties of the covariance kernel */
template <class Space>
struct kernel_traits<covariance_kernel<Space> >
{
	typedef typename build<location<Space> >::type test_input_t;
	typedef typename build<location<Space> >::type trial_input_t;
	typedef collect<covariance_data> data_t;
	typedef typename single_brick_wall<CKernel<Space> >::type output_t;
	enum { result_rows = 1, result_cols = 1 };
	typedef gauss_family_tag quadrature_family_t;
	static bool const is_symmetric = true;
	typedef asymptotic::inverse<1> far_field_behaviour_t;
	static bool const is_singular = false;
};

template <class Space>
class covariance_kernel :
	public kernel_base<covariance_kernel<Space> >
{
public:
	covariance_kernel(double sigma, double length) :
		kernel_base<covariance_kernel>(covariance_data(sigma, length)) {}
};

}

#endif // COVARIANCE_KERNEL_HPP_INCLUDED
