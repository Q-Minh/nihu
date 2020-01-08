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
 *
 * \todo Check non-stationary implementation options
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
	typedef Eigen::Matrix<double, field_dimension, field_dimension> field_variance_t;
	
	exponential_covariance_kernel(field_variance_t const &variance, double length) :
		m_field_variance(variance), m_length(length)
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
	
	field_variance_t const &get_variance(void) const
	{
		return m_field_variance;
	}
	
	double get_correlation_length(void) const
	{
		return m_length;
	}

private:
	field_variance_t m_field_variance;	/**< \brief Variance matrix */
	double m_length;		/**< \brief Correlation length */
};



template <class Space, class Dimension>
class gaussian_covariance_kernel;

/// Gaussian covariance kernel traits
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

template <class Space>
class ConstVariance
{
public:
	typedef Space space_t;
	typedef typename space_t::scalar_t scalar_t;
	static size_t const space_dimension = Space::dimension;
	typedef typename space_t::location_t location_t;

	typedef Eigen::Matrix<scalar_t, space_dimension, space_dimension> result_t;

	ConstVariance(result_t const &variance)
		: m_variance(variance)
	{
	}

	result_t const &get_variance(location_t const&) const
	{
		return m_variance;
	}

private:
	result_t m_variance;
};

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
	typedef Eigen::Matrix<double, space_dimension, space_dimension> space_variance_t;
	typedef Eigen::Matrix<double, field_dimension, field_dimension> field_variance_t;

	gaussian_covariance_kernel(
		field_variance_t const &field_variance, 
		space_variance_t const &space_variance
	) 
		: m_field_variance(field_variance)
		, m_space_variance(space_variance)
		, m_inv_space_variance(m_space_variance.inverse())
	{
	}

	result_t operator()(location_t const &x, location_t const &y) const
	{
		location_t d = x - y;
		double Q = d.transpose() * m_inv_space_variance * d;
		return get_field_variance() * std::exp(-Q);
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}

	field_variance_t const &get_field_variance(void) const
	{
		return m_field_variance;
	}

	space_variance_t const &get_space_variance(void) const
	{
		return m_space_variance;
	}

private:
	field_variance_t m_field_variance;		/**< \brief Field variance matrix */
	space_variance_t m_space_variance;		/**< \brief Space variance matrix */
	space_variance_t m_inv_space_variance;	/**< \brief inverse of the space variance matrix */
};

#if 0
// NOTE: NON-STATIONARY STARTED HERE

template <class SpaceVariance, class FieldVariance>
class gaussian_covariance_kernel;

/// Gaussian covariance kernel traits
namespace kernel_traits_ns
{
template <class SpaceVariance, class FieldVariance>
struct space<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > : typename SpaceVariance::space_t {};

template<class SpaceVariance, class FieldVariance>
struct result<gaussian_covariance_kernel<SpaceVariance, FieldVariance> >
{
	typedef FieldVariance::result_t type;
};

template <class SpaceVariance, class FieldVariance>
struct quadrature_family<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > : gauss_family_tag {};

template <class SpaceVariance, class FieldVariance>
struct is_singular<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > : std::false_type {};

template <class SpaceVariance, class FieldVariance>
struct test_input<gaussian_covariance_kernel<SpaceVariance, FieldVariance> >
{
	typedef location_input<typename SpaceVariance::space_t> type;
};

template <class SpaceVariance, class FieldVariance>
struct trial_input<gaussian_covariance_kernel<SpaceVariance, FieldVariance> >
{
	typedef location_input<typename SpaceVariance::space_t> type;
};

template <class SpaceVariance, class FieldVariance>
struct is_symmetric<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > : std::true_type {};

/** \todo this is incorrect but works :) */
template <class SpaceVariance, class FieldVariance>
struct far_field_behaviour<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > : asymptotic::inverse<1> {};
} // end of NiHu::kernel_traits_ns

template <class SpaceVariance, class FieldVariance>
class gaussian_covariance_kernel :
	public kernel_base<gaussian_covariance_kernel<SpaceVariance, FieldVariance> >
{
	typedef SpaceVariance space_variance_t;
	typedef FieldVariance field_variance_t;

	typedef kernel_base<gaussian_covariance_kernel<SpaceVariance, FieldVariance> > base_t;
	typedef typename base_t::test_input_t test_input_t;
	typedef typename base_t::trial_input_t trial_input_t;

	typedef typename field_variance_t::result_t field_var_t;
	typedef typename space_variance_t::result_t space_var_t;

	gaussian_covariance_kernel(field_variance_t const &fvar, space_variance_t const &svar)
		: m_field_variance(fvar)
		, m_space_variance(svar)
	{
	}

	result_t operator()(location_t const &x, location_t const &y) const
	{
		location_t d = x - y;
		space_variance_t svar_x = get_space_variance(x);
		space_variance_t svar_y = get_space_variance(y);
		field_variance_t fvar = get_field_variance(x);

		space_variance_t svar = .5 * (svar_x + svar_y);
		double Q = d.transpose() * svar.inverse() * d;
		double scale = std::pow(svar_x.determinant(), .25) * std::pow(svar_y.determinant(), .25) / std::sqrt(svar.determinant());

		return fvar * scale * std::exp(-Q);
	}

	field_var_t get_field_variance(location_t const &x) const
	{
		return m_field_variance.get_variance(x);
	}

	space_var_t get_space_variance(location_t const &x) const
	{
		return m_space_variance.get_variance(x);
	}

	result_t operator()(
		test_input_t const &x,
		trial_input_t const &y) const
	{
		return (*this)(x.get_x(), y.get_x());
	}

private:
	space_variance_t m_space_variance;
	field_variance_t m_field_variance;
};
// NOTE: NON-STATIONARY ENDS HERE
#endif

} // end of namespace NiHu

#endif // COVARIANCE_KERNEL_HPP_INCLUDED
