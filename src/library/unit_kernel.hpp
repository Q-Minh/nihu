// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
 * \file unit_kernel.hpp
 * \ingroup library
 * \brief implementation of the unit kernel \f$K(x,y) = 1\f$
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef UNIT_KERNEL_HPP_INCLUDED
#define UNIT_KERNEL_HPP_INCLUDED

#include "../core/kernel.hpp"
#include "../util/brick.hpp"
#include "../core/gaussian_quadrature.hpp"

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
template <class Scalar>
class unit_kernel;

/** \brief traits of the unit kernel */
template<class Scalar>
struct kernel_traits<unit_kernel<Scalar> >
{
	/** \brief test input type */
	typedef empty_wall test_input_t;
	/** \brief trial input type */
	typedef empty_wall trial_input_t;
	/** \brief kernel result type */
	typedef Scalar result_t;
	/** \brief quadrature family tag */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief shows if kernel is symmetric */
	static bool const is_symmetric = true;
	/** \brief the singularity order (r^-order) */
	static unsigned const singularity_order = 0;
	/** \brief the singular quadrature order */
	static unsigned const singular_quadrature_order = 0;
};

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
template <class Scalar>
class unit_kernel :
	public kernel_base<unit_kernel<Scalar> >
{
public:
	/** \brief the crtp base type */
	typedef kernel_base<unit_kernel<Scalar> > base_t;
	/** \brief the scalar type */
	typedef typename base_t::scalar_t scalar_t;
	/** \brief the result type */
	typedef typename base_t::result_t result_t;
	/** \brief the test input type */
	typedef typename base_t::test_input_t test_input_t;
	/** \brief the trial input type */
	typedef typename base_t::trial_input_t trial_input_t;

	/**
	 * \brief evaluate kernel at test and trial positions
	 * \return kernel value K(x,y)
	 */
	constexpr result_t operator()(
		test_input_t const &, trial_input_t const &) const
	{
		return result_t(1.0);
	}

	/**
	 * \brief evaluate kernel complexity
	 * \return kernel value K(x,y)
	 */
	constexpr unsigned estimate_complexity(
		test_input_t const &, trial_input_t const &, scalar_t const &) const
	{
		return 0;
	}
};

#endif // UNIT_KERNEL_HPP_INCLUDED

