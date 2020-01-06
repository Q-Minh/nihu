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
 * \brief Implementation of kernels representing covariance functions of stochastic processes
 * \author Peter Fiala fiala@hit.bme.hu 
 * \author Peter Rucz rucz@hit.bme.hu
 * \ingroup lib_kernel
 */

#ifndef COVARIANCE_KERNEL_HPP_INCLUDED
#define COVARIANCE_KERNEL_HPP_INCLUDED

#include "../core/field_dimension.hpp"
#include "../core/global_definitions.hpp"
#include "../core/kernel.hpp"
#include "../core/gaussian_quadrature.hpp"
#include "location_normal.hpp"

namespace NiHu
{

template <class Space, class Dimension>
class exponential_covariance_kernel;

/// GENERAL TRAITS
namespace kernel_traits_ns
{
	template <class Space, class Dimension>
	struct space<exponential_covariance_kernel<Space, Dimension> > : Space {};

	template<class Space, class Dimension>
	struct result<exponential_covariance_kernel<Space, Dimension> >
	{
		typedef Eigen::Matrix<typename Space::scalar_t, Dimension::value, Dimension::value> type;
	};
	
	// Specialisation
	template<class Space>
	struct result<exponential_covariance_kernel<Space, field_dimension::_1d> >
	{
		typedef typename Space::scalar_t type;
	};

	template <class Space, class Dimension>
	struct quadrature_family<exponential_covariance_kernel<Space, Dimension> > : gauss_family_tag {};

	template <class Space, class Dimension>
	struct is_singular<exponential_covariance_kernel<Space, Dimension> > : std::false_type {};

	template <class Space, class Dimension>
	struct test_input<exponential_covariance_kernel<Space, Dimension> >
	{
		typedef location_input<Space> type;
	};

	template <class Space, class Dimension>
	struct trial_input<exponential_covariance_kernel<Space, Dimension> >
	{
		typedef location_input<Space> type;
	};

	template <class Space, class Dimension>
	struct is_symmetric<exponential_covariance_kernel<Space, Dimension> > : std::true_type {};

	/** \todo this is incorrect but works :) */
	template <class Space, class Dimension>
	struct far_field_behaviour<exponential_covariance_kernel<Space, Dimension> > : asymptotic::inverse<1> {};
}


template <class Space, class Dimension>
class exponential_covariance_kernel :
	public kernel_base<exponential_covariance_kernel<Space, Dimension> >
{
public:
	typedef kernel_base<exponential_covariance_kernel<Space, Dimension> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::scalar_t distance_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t location_t;
	typedef Dimension dimension_t;
	size_t static const field_dimension = dimension_t::value;
	
	exponential_covariance_kernel(result_t const &variance, double length) :
		m_variance(variance), m_length(length)
	{
	}
	
	static distance_t distance(location_t const &x, location_t const &y)
	{
		return (x-y).norm();
	}
	
	static distance_t distance_sph(location_t const &x, location_t const &y)
	{
		result_t v = x.normalized().dot(y.normalized());
		v = std::min(v, 1.0);
		v = std::max(v, -1.0);
		return std::acos(v);
	}
	
	result_t operator()(location_t const &x, location_t const &y) const
	{
		auto r = distance(x, y);
		return get_variance() * std::exp(-r / get_correlation_length());
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}
	
	result_t const &get_variance(void) const
	{
		return m_variance;
	}
	
	double get_correlation_length(void) const
	{
		return m_length;
	}

private:
	result_t m_variance;	/**< \brief Variance matrix */
	double m_length;		/**< \brief Correlation length */
};



template <class Space, class Dimension>
class gaussian_covariance_kernel;

/// GENERAL TRAITS
namespace kernel_traits_ns
{
template <class Space, class Dimension>
struct space<gaussian_covariance_kernel<Space, Dimension> > : Space {};

template<class Space, class Dimension>
struct result<gaussian_covariance_kernel<Space, Dimension> >
{
	typedef Eigen::Matrix<typename Space::scalar_t, Dimension::value, Dimension::value> type;
};

// Specialisation
template<class Space>
struct result<gaussian_covariance_kernel<Space, field_dimension::_1d> >
{
	typedef typename Space::scalar_t type;
};

template <class Space, class Dimension>
struct quadrature_family<gaussian_covariance_kernel<Space, Dimension> > : gauss_family_tag {};

template <class Space, class Dimension>
struct is_singular<gaussian_covariance_kernel<Space, Dimension> > : std::false_type {};

template <class Space, class Dimension>
struct test_input<gaussian_covariance_kernel<Space, Dimension> >
{
	typedef location_input<Space> type;
};

template <class Space, class Dimension>
struct trial_input<gaussian_covariance_kernel<Space, Dimension> >
{
	typedef location_input<Space> type;
};

template <class Space, class Dimension>
struct is_symmetric<gaussian_covariance_kernel<Space, Dimension> > : std::true_type {};

/** \todo this is incorrect but works :) */
template <class Space, class Dimension>
struct far_field_behaviour<gaussian_covariance_kernel<Space, Dimension> > : asymptotic::inverse<1> {};
}


template <class Space, class Dimension>
class gaussian_covariance_kernel :
	public kernel_base<gaussian_covariance_kernel<Space, Dimension> >
{
public:
	typedef kernel_base<gaussian_covariance_kernel<Space, Dimension> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::scalar_t distance_t;
	typedef typename base_t::result_t result_t;
	typedef typename base_t::x_t location_t;
	typedef Dimension dimension_t;
	size_t static const field_dimension = dimension_t::value;
	size_t static const space_dimension = Space::dimension;
	typedef Eigen::Matrix<double, space_dimension, space_dimension> sigma_t;

	gaussian_covariance_kernel(result_t const &variance, sigma_t const &sigma) :
		m_variance(variance), m_sigma(sigma), m_inv_sigma(m_sigma.inverse())
	{
	}

	result_t operator()(location_t const &x, location_t const &y) const
	{
		location_t d = x - y;
		double Q = d.transpose() * m_inv_sigma * d;
		return get_variance() * std::exp(-Q);
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}

	result_t const &get_variance(void) const
	{
		return m_variance;
	}

	sigma_t const &get_sigma(void) const
	{
		return m_sigma;
	}

private:
	result_t m_variance;	/**< \brief Variance matrix */
	sigma_t m_sigma;		/**< \brief Variance matrix of the correlation */
	sigma_t m_inv_sigma;	/**< \brief inverse of the variance matrix */
};


} // end of namespace NiHu

#endif // COVARIANCE_KERNEL_HPP_INCLUDED
