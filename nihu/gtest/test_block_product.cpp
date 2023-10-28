#include <gtest/gtest.h>
#include "nihu/util/block_product.hpp"

#include <type_traits>

TEST(BlockProduct, EigenDouble)
{
	typedef Eigen::Matrix<double, 2, 1> left_t;
	typedef Eigen::Matrix<double, 3, 3> mat_t;
	typedef Eigen::Matrix<double, 4, 1> right_t;

	typedef NiHu::block_product_result_type<left_t, mat_t, right_t>::type brt_t;

	static_assert(std::is_same<brt_t, Eigen::Matrix<double, 2*3, 3*4> >::value,
		"Block product with EigenDouble failed");

	left_t left(left_t::Constant(1.0));
	mat_t mat(mat_t::Constant(2.0));
	right_t right(right_t::Constant(3.0));

	typedef decltype(mat_t::Constant(1.0)) tt;

	bool const b = std::is_base_of<Eigen::MatrixBase<std::decay<tt>::type>, std::decay<tt>::type>::value;
	std::cout << b << std::endl;

	bool const c = std::is_base_of<Eigen::MatrixBase<tt>, tt>::value;
	std::cout << c << std::endl;


	brt_t brt = NiHu::block_product(left_t::Constant(1.0), mat_t::Constant(2.0), right_t::Constant(3.0));

	std::cout << brt << std::endl;
}

#if 0
TEST(BlockProduct, EigenComplex)
{
	typedef Eigen::Matrix<double, 1, 1> left;
	typedef Eigen::Matrix<std::complex<double>, 1, 1> mat;
	typedef Eigen::Matrix<double, 1, 1> right;

	typedef NiHu::block_product_result_type<left, mat, right>::type brt1_t;

	brt1_t brt1 = NiHu::block_product(
		left::Constant(1.0), mat::Constant(2.0), right::Constant(3.0));
}

TEST(BlockProduct, Double)
{
	typedef Eigen::Matrix<double, 1, 1> left;
	typedef double mat;
	typedef Eigen::Matrix<double, 1, 1> right;

	typedef NiHu::block_product_result_type<left, mat, right>::type brt1_t;

	brt1_t brt1 = NiHu::block_product(
		left::Constant(1.0), 2.0, right::Constant(3.0));
}

TEST(BlockProduct, Complex)
{
	typedef Eigen::Matrix<double, 1, 1> left;
	typedef std::complex<double> mat;
	typedef Eigen::Matrix<double, 1, 1> right;

	typedef NiHu::block_product_result_type<left, mat, right>::type brt1_t;

	brt1_t brt1 = NiHu::block_product(
		left::Constant(1.0), 2.0, right::Constant(3.0));
}
#endif