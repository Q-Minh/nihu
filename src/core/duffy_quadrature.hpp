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
 * \file duffy_quadrature.hpp
 * \ingroup quadrature
 * \brief Duffy-type singular quadrature transformations
 */

#ifndef DUFFY_QUADRATURE_HPP_INCLUDED
#define DUFFY_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

/** \brief Traits class of a Duffy quadrature */
template <class LSet>
struct duffy_traits;

/**
 * \brief Transform regular quadratures into weakly singular Duffy-type quadratures
 * \tparam QuadFamily the regular quadrature family
 * \tparam LSet the element geometrical shape set representation
 */
template <class QuadFamily, class LSet>
class duffy_quadrature
{
	// CRTP check
	static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
		"The LSet parameter must be derived from shape_set_base<LSet>");
public:
	/** \brief template argument as nested type */
	typedef QuadFamily quadrature_family_tag;
	/** \brief template argument as nested type */
	typedef LSet lset_t;

	/** \brief element's domain type */
	typedef typename lset_t::domain_t domain_t;
	/** \brief base domain variable type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief domain's scalar type */
	typedef typename domain_t::scalar_t scalar_t;
	/** \brief the return quadrature type */
	typedef typename quadrature_type<quadrature_family_tag, domain_t>::type ret_quad_type;

private:
	/** \brief the source quadrature type (regular quad quadrature)
		\todo quad domain is hard coded here
	*/
	typedef typename quadrature_type<quadrature_family_tag, quad_domain>::type source_quad_type;
	/** \brief the local coordinate matrix type of the transformation */
	typedef Eigen::Matrix<scalar_t, quad_domain::num_corners, domain_t::dimension> coords_t;

	/**
	 * \brief Helper function for Duffy quadrature generation
	 * \param [in] degree the polynomial degree of the original regular quadrature
	 * \param [in] array indices of the Duffy triangle's nonsingular corners
	 * \param [in] singular_point coordinates of the singular point
	 * \param [out] reference to the result quadrature
	 * \return reference to the result for chaining
	 * \todo quad_1_shape_set is hardcoded
	 */
	static ret_quad_type &duffy_impl(
		unsigned degree,
		unsigned const *array,
		xi_t const &sing_coord,
		ret_quad_type &result)
	{
		unsigned num_duffies = *array;
		array++;

		source_quad_type source(degree);

		coords_t coords;
		coords.row(0) = coords.row(1) = sing_coord;

		for (size_t d = 0; d < num_duffies; ++d)
		{
			coords.row(2) = domain_t::get_corner(array[d]);
			coords.row(3) = domain_t::get_corner(array[d+1]);

			result += source.template transform<quad_1_shape_set>(coords);
		}

		return result;
	}

public:
	/**
	 * \brief create a duffy quadrature that is singular on one corner of the selected element
	 * \param [in] degree the polynomial degree of the original regular quadrature
	 * \param [in] singular_corner index of the singular corner
	 * \return a Duffy type singular quadrature
	 */
	static ret_quad_type on_corner(unsigned degree, unsigned singular_corner)
	{
		ret_quad_type result;
		return duffy_impl(
			degree,
			duffy_traits<lset_t>::duffy_corner_indices[singular_corner],
			lset_t::corner_at(singular_corner),
			result);
	}

	/**
	 * \brief create a duffy quadrature that is singular on the face of the selected element
	 * \param [in] degree the polynomial degree of the original regular quadrature
	 * \param [in] singular_point coordinates of the singular point
	 * \return a Duffy type singular quadrature
	 */
	static ret_quad_type on_face(unsigned degree, xi_t const &singular_point)
	{
		static unsigned const n = domain_t::num_corners;
		unsigned face_indices[n+2];
		face_indices[0] = n;
		for (unsigned c = 0; c <= n; ++c)
			face_indices[1+c] = c % n;
		ret_quad_type result;
		return duffy_impl(degree, face_indices, singular_point, result);
	}
}; // end of class duffy_quadrature



/** \brief Specialisation of ::duffy_traits for ::tria_1_shape_set */
template <>
struct duffy_traits<tria_1_shape_set>
{
	/** \brief indices of the Duffy corners for singular corners */
	static unsigned const duffy_corner_indices[3][2+1];
};

unsigned const duffy_traits<tria_1_shape_set>::duffy_corner_indices[3][2+1] = {
	{1, /*|*/ 1, 2},
	{1, /*|*/ 2, 0},
	{1, /*|*/ 0, 1}
};


/** \brief Specialisation of ::duffy_traits for ::quad_1_shape_set */
template <>
struct duffy_traits<quad_1_shape_set>
{
	/** \brief indices of the Duffy corners for singular corners */
	static unsigned const duffy_corner_indices[4][3+1];
	/** \brief indices of the Duffy corners for internal singular point */
	static unsigned const duffy_face_indices[5+1];
};

unsigned const duffy_traits<quad_1_shape_set>::duffy_corner_indices[4][3+1] = {
	{2, /*|*/ 1, 2, 3},
	{2, /*|*/ 2, 3, 0},
	{2, /*|*/ 3, 0, 1},
	{2, /*|*/ 0, 1, 2}
};


/** \brief Specialisation of ::duffy_traits for ::tria_2_shape_set */
template <>
struct duffy_traits<tria_2_shape_set>
{
	/** \brief indices of the Duffy corners for singular corners */
	static unsigned const duffy_corner_indices[6][3+1];
};

unsigned const duffy_traits<tria_2_shape_set>::duffy_corner_indices[6][3+1] = {
	{1, /*|*/ 1, 2},
	{2, /*|*/ 1, 2, 0},
	{1, /*|*/ 2, 0},
	{2, /*|*/ 2, 0, 1},
	{1, /*|*/ 0, 1},
	{2, /*|*/ 0, 1, 2}
};


/** \brief Specialisation of ::duffy_traits for ::quad_2_shape_set */
template <>
struct duffy_traits<quad_2_shape_set>
{
	/** \brief indices of the Duffy corners for singular corners */
	static unsigned const duffy_corner_indices[9][5+1];
};

unsigned const duffy_traits<quad_2_shape_set>::duffy_corner_indices[9][5+1] = {
	{2, /*|*/ 1, 2, 3},
	{3, /*|*/ 1, 2, 3, 0},
	{2, /*|*/ 2, 3, 0},
	{3, /*|*/ 2, 3, 0, 1},
	{2, /*|*/ 3, 0, 1},
	{3, /*|*/ 3, 0, 1, 2},
	{2, /*|*/ 0, 1, 2},
	{3, /*|*/ 0, 1, 2, 3},
	{4, /*|*/ 0, 1, 2, 3, 0}
};

#endif // DUFFY_QUADRATURE_HPP_INCLUDED

