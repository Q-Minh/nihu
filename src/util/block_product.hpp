#ifndef BLOCK_PRODUCT_HPP_INCLUDED
#define BLOCK_PRODUCT_HPP_INCLUDED

#include "eigen_utils.hpp"

namespace internal
{
	template <class left, class mat, class right>
	class block_product_impl;

	template <class scalar1, class scalar2, class scalar3, int N1, int N2, int N3>
	class block_product_impl<
		Eigen::Matrix<scalar1, N1, 1>,
		Eigen::Matrix<scalar2, N2, N2>,
		Eigen::Matrix<scalar3, 1, N3>
	>
	{
	public:
		typedef Eigen::Matrix<scalar1, N1*N2, N2*N3> result_type;

		static result_type eval(
			Eigen::Matrix<scalar1, N1, 1> const &v1,
			Eigen::Matrix<scalar2, N2, N2> const &m,
			Eigen::Matrix<scalar3, 1, N3> const &v2)
		{
			result_type result;
			for (int row = 0; row < N1; ++row)
				for (int col = 0; col < N3; ++col)
					result.block(row*N2, col*N2, N2, N2) = v1(row) * m * v2(col);
			return result;
		}
	};
}


template <class left, class mat, class right>
typename internal::block_product_impl<left, mat, right>::result_type
block_product(left const &l, mat const &m, right const &r)
{
	return internal::block_product_impl<left, mat, right>::eval(l, m, r);
}

#endif // BLOCK_PRODUCT_HPP_INCLUDED
