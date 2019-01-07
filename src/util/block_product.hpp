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
 * \brief declaration of template function ::block_product, ::semi_block_product and its related metafunctions
 */

#ifndef BLOCK_PRODUCT_HPP_INCLUDED
#define BLOCK_PRODUCT_HPP_INCLUDED

#include "eigen_utils.hpp"
#include "plain_type.hpp"
#include "product_type.hpp"

#include <iostream>

namespace NiHu
{

namespace internal
{
template <class leftDerived, class mat, class rightDerived, 
	bool isEigen = is_eigen<mat>::value>
class block_product_impl;

/** \brief specialisation of block_product for two Eigen vectors and something general */
template <class leftDerived, class Scalar, class rightDerived>
class block_product_impl<leftDerived, Scalar, rightDerived, false>
{
public:
	typedef typename plain_type<
		typename product_type<
		Scalar,
		typename plain_type<
		typename product_type<leftDerived, Eigen::Transpose<rightDerived> >::type
		>::type
		>::type
	>::type result_type;

	static result_type eval(
		Eigen::MatrixBase<leftDerived> const &v1,
		Scalar const &m,
		Eigen::MatrixBase<rightDerived> const &v2
	)
	{
		return m * (v1 * v2.transpose());
	}
};

/** \brief specialisation of block_product for three Eigen matrices */
template <class leftDerived, class matDerived, class rightDerived>
class block_product_impl<leftDerived, matDerived, rightDerived, true>
{
	typedef typename leftDerived::Scalar scalar1;
	typedef typename matDerived::Scalar scalar2;
	typedef typename rightDerived::Scalar scalar3;
	enum {
		N1 = leftDerived::RowsAtCompileTime,
		N21 = matDerived::RowsAtCompileTime, N22 = matDerived::ColsAtCompileTime,
		N3 = rightDerived::RowsAtCompileTime
	};
	typedef typename product_type<
		scalar1,
		typename product_type<scalar2, scalar3>::type
	>::type scalar;

public:
	typedef Eigen::Matrix<scalar, N1*N21, N22*N3> result_type;

	static result_type eval(
		Eigen::MatrixBase<leftDerived> const &v1,
		Eigen::MatrixBase<matDerived> const &m,
		Eigen::MatrixBase<rightDerived> const &v2)
	{
		result_type result;
		for (size_t row = 0; row < N1; ++row)
			for (size_t col = 0; col < N3; ++col)
				result.template block<N21, N22>(row*N21, col*N22) =
				m * static_cast<scalar>(v1(row, 0) * v2(col, 0));
		return result;
	}
};


template <class mat, class right>
class semi_block_product_impl;

template <class mat, class right>
class semi_block_product_impl
{
public:
	typedef typename plain_type<
		typename product_type<
		mat,
		Eigen::Transpose<right>
		>::type
	>::type result_type;

	static result_type eval(mat const &m, Eigen::MatrixBase<right> const &v2)
	{
		return m * v2.transpose();
	}
};


template <class scalar2, int N2, class right>
class semi_block_product_impl<Eigen::Matrix<scalar2, N2, N2>, right>
{
	typedef typename right::Scalar scalar3;
	enum { N3 = right::RowsAtCompileTime };
	typedef typename product_type<scalar2, scalar3>::type scalar;
public:
	typedef Eigen::Matrix<scalar, N2, N2*N3> result_type;

	static result_type eval(
		Eigen::MatrixBase<Eigen::Matrix<scalar2, N2, N2> > const &m,
		Eigen::MatrixBase<right> const &v2)
	{
		result_type result;
		for (int col = 0; col < N3; ++col)
			result.template block<N2, N2>(0, col*N2) = m * v2(col);
		return result;
	}
};

} // end of namespace internal



/** \brief metafunction returning the value type of a block product */
template <class leftDerived, class mat, class rightDerived>
struct block_product_result_type
{
	typedef typename internal::block_product_impl<leftDerived, mat, rightDerived>::result_type type;
};

/** \brief compute a block product l * m * r^T
* \tparam leftDerived the left Eigen vector type
* \tparam mat the matrix type
* \tparam rightDerived the right Eigen vector type
* \param [in] l the left Eigen vector
* \param [in] m the matrix
* \param [in] r the right Eigen vector
* \return the block product l * m * r^T
*/
template <class leftDerived, class mat, class rightDerived>
typename block_product_result_type<leftDerived, mat, rightDerived>::type
block_product(
	Eigen::MatrixBase<leftDerived> const &l,
	mat const &m,
	Eigen::MatrixBase<rightDerived> const &r)
{
	return internal::block_product_impl<leftDerived, mat, rightDerived>::eval(l, m, r);
}

/** \brief compute semi block product of a matrix and a vector m * v^T
 * \tparam mat the matrix type
 * \tparam right the Eigen vector type
 * \param [in] m the matrix
 * \param [in] r the Eigen vector
 * \return the block product m * r^T
 */
template <class mat, class right>
auto semi_block_product(mat const &m, Eigen::MatrixBase<right> const &r)
-> decltype(internal::semi_block_product_impl<mat, right>::eval(m, r))
{
	return internal::semi_block_product_impl<mat, right>::eval(m, r);
}

/** \brief metafunction returning the value type of a semi block product */
template <class mat, class right>
struct semi_block_product_result_type
{
	typedef typename internal::semi_block_product_impl<mat, right>::result_type type;
};

} // end of namespace NiHu

#endif // BLOCK_PRODUCT_HPP_INCLUDED

