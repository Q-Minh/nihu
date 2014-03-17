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

/** \file block_product.hpp
 * \ingroup util
 * \brief declaration of template function ::block_product and its related metafunctions
 */

#ifndef BLOCK_PRODUCT_HPP_INCLUDED
#define BLOCK_PRODUCT_HPP_INCLUDED

#include "eigen_utils.hpp"
#include "product_type.hpp"
#include "plain_type.hpp"
#include <iostream>

namespace internal
{
	template <class left, class mat, class right>
	class block_product_impl;

	/** \brief specialisation of block_product for two Eigen vectors and something general */
	template <class scalar1, class mat, class scalar3, int N1, int N3>
	class block_product_impl<
		Eigen::Matrix<scalar1, N1, 1>,
		mat,
		Eigen::Matrix<scalar3, N3, 1>
	>
	{
		typedef Eigen::Matrix<scalar1, N1, 1> lhs_t;
		typedef Eigen::Matrix<scalar3, N3, 1> rhs_t;

	public:
		typedef typename plain_type<
			typename product_type<
				mat,
				typename plain_type<
					typename product_type<
						lhs_t, Eigen::Transpose<rhs_t>
					>::type
				>::type
			>::type
		>::type result_type;

		static result_type eval(
			Eigen::MatrixBase<lhs_t> const &v1,
			mat const &m,
			Eigen::MatrixBase<rhs_t> const &v2
		)
		{
			return m * (v1 * v2.transpose());
		}
	};


	/** \brief specialisation of block_product for three Eigen matrices */
	template <class scalar1, class scalar2, class scalar3, int N1, int N2, int N3>
	class block_product_impl<
		Eigen::Matrix<scalar1, N1, 1>,
		Eigen::Matrix<scalar2, N2, N2>,
		Eigen::Matrix<scalar3, N3, 1>
	>
	{
		typedef typename product_type<
			scalar1,
			typename product_type<scalar2, scalar3>::type
		>::type scalar;

	public:
		typedef Eigen::Matrix<scalar, N1*N2, N2*N3> result_type;

		static result_type eval(
			Eigen::MatrixBase<Eigen::Matrix<scalar1, N1, 1> > const &v1,
			Eigen::MatrixBase<Eigen::Matrix<scalar2, N2, N2> > const &m,
			Eigen::MatrixBase<Eigen::Matrix<scalar3, N3, 1> > const &v2)
		{
			result_type result;
			for (Index row = 0; row < N1; ++row)
				for (Index col = 0; col < N3; ++col)
					result.block(row*N2, col*N2, N2, N2) = v1(row) * m * v2(col);
			return result;
		}
	};
}

/** \brief evaluate the block product v1 * m * v2^T
 * \tparam v1 the left hand side vector type
 * \tparam mat the matrix type
 * \tparam v2 the right hand side vector type
 * \param l the left vector
 * \param m the center matrix
 * \param r the right vector
 */
template <class v1, class mat, class v2>
auto block_product(Eigen::MatrixBase<v1> const &l, mat const &m, Eigen::MatrixBase<v2> const &r) ->
	decltype(internal::block_product_impl<v1, mat, v2>::eval(l, m, r))
{
	return internal::block_product_impl<v1, mat, v2>::eval(l, m, r);
}

/** \brief metafunction returning the value type of a block product */
template <class v1, class mat, class v2>
struct block_product_result_type
{
	typedef typename internal::block_product_impl<v1, mat, v2>::result_type type;
};

#endif // BLOCK_PRODUCT_HPP_INCLUDED

