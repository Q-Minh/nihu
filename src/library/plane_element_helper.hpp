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

/** \file plane_element_helper.hpp
 * \brief helper functions to compute analytical integrals over plane elements
 */

#ifndef PLANE_ELEMENT_HELPER_HPP_INCLUDED
#define PLANE_ELEMENT_HELPER_HPP_INCLUDED

namespace NiHu
{

/**
 * \brief compute angles and radii in a plane element
 * \param [in] elem the element
 * \param [out] r the radii from the collocation point
 * \param [out] theta the angles between the radii and the sides
 * \param [out] alpha the angles between the radii
 */
template <class elem_t>
void plane_element_helper(
	elem_t const &elem,
	typename elem_t::x_t const &x0,
	typename elem_t::scalar_t r[],
	typename elem_t::scalar_t theta[],
	typename elem_t::scalar_t alpha[])
{
	enum{ N = elem_t::domain_t::num_corners };

	auto const &C_old = elem.get_coords();

	typename elem_t::coords_t R, C;
	for (unsigned i = 0; i < N; ++i)
	{
		R.col(i) = C_old.col(i) - x0;
		r[i] = R.col(i).norm();
		R.col(i) /= r[i];
		C.col(i) = C_old.col(i) - C_old.col((i+1) % N);
		C.col(i) /= C.col(i).norm();
	}

	for (unsigned i = 0; i < N; ++i)
	{
		theta[i] = std::acos(R.col(i).dot(R.col((i+1) % N)));
		alpha[i] = std::acos(R.col(i).dot(C.col(i)));
	}
}

}

#endif // PLANE_ELEMENT_HELPER_HPP_INCLUDED
