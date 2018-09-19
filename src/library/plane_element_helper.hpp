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
 * \param [out] theta the angles between the radii
 * \param [out] alpha the angles between the radii and the sides
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
		double cs = R.col(i).dot(R.col((i+1) % N));
		if (cs >= 1.)
			theta[i] = 0.;
		else
			theta[i] = std::acos(cs);
		alpha[i] = std::acos(R.col(i).dot(C.col(i)));
	}
}

template <class V>
Eigen::Matrix<double, 3, 3> plane_elem_transform(
		Eigen::DenseBase<V> const &v1_in, 
		Eigen::DenseBase<V> const &v2_in)
{
	// \todo optimize out these instantiations please
	Eigen::Matrix<double, 3, 1> v1 = v1_in;
	Eigen::Matrix<double, 3, 1> v2 = v2_in;
	
	Eigen::Matrix<double, 3, 3> T;
	T.col(0) = v1.normalized();
	T.col(2) = v1.cross(v2).normalized();
	T.col(1) = T.col(2).cross(T.col(0)).normalized();
	
	return T;
}

} // end of namespace NiHu

#endif // PLANE_ELEMENT_HELPER_HPP_INCLUDED
