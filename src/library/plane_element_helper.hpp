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

/** 
 * \file plane_element_helper.hpp
 * \brief helper functions to compute analytical integrals over plane elements
 * \ingroup library
 */

#ifndef NIHU_PLANE_ELEMENT_HELPER_HPP_INCLUDED
#define NIHU_PLANE_ELEMENT_HELPER_HPP_INCLUDED

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
	/** \todo check if elem_t::num_corners here */
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

template <class matrixDerived, class vectorDerived>
void plane_elem_helper_mid(
	Eigen::DenseBase<matrixDerived> const &coords,
	Eigen::DenseBase<vectorDerived> const &x0,
	double ref_distance[],
	double theta_lim[],
	double theta0[]
)
{
	// geometrical parameters (planar helpers)
	unsigned const N = coords.cols();
	for (unsigned n = 0; n < N; ++n)
	{
		Eigen::Matrix<double, 2, 1> c1 = coords.topRows(2).col(n);			// corner
		Eigen::Matrix<double, 2, 1> c2 = coords.topRows(2).col((n + 1) % N);	// next corner
		
		Eigen::Matrix<double, 2, 1> l = (c2 - c1).normalized();					// side unit vector

		Eigen::Matrix<double, 2, 1> d1 = c1 - x0.topRows(2);			// vector to corners
		Eigen::Matrix<double, 2, 1> d0 = d1 - l * d1.dot(l);		// perpendicular to side

		theta_lim[n] = std::atan2(d1(1), d1(0));	// corner angle
		theta0[n] = std::atan2(d0(1), d0(0));	// mid angle
		ref_distance[n] = d0.norm();			// distance to side
	}
}

/**
 * \brief Transformation matrix to get planar element parallel to the x-y plane
 * \tparam V Vector type, must be an Eigen vector
 * \returns Transformation matrix of 3x3 type
 * \details
 * Using the transformation matrix, the element can be projected such that it 
 * becomes parallel with the x-y plane. Note that as there is no translation in 
 * the transform, the element will not necessarily be shifted into the z=0 
 * plane.
 */
template <class V>
Eigen::Matrix<double, 3, 3> plane_elem_transform(
		Eigen::DenseBase<V> const &v1_in, 
		Eigen::DenseBase<V> const &v2_in)
{
	/** \todo optimize out these instantiations please */
	Eigen::Matrix<double, 3, 1> v1 = v1_in;
	Eigen::Matrix<double, 3, 1> v2 = v2_in;
	
	Eigen::Matrix<double, 3, 3> T;
	T.col(0) = v1.normalized();
	T.col(2) = v1.cross(v2).normalized();
	T.col(1) = T.col(2).cross(T.col(0)).normalized();
	
	return T;
}

} // end of namespace NiHu

#endif /* NIHU_PLANE_ELEMENT_HELPER_HPP_INCLUDED */
