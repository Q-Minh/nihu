#include <gtest/gtest.h>
#include "library/lib_shape.hpp"

// testing line_0_shape_set
static_assert(std::is_same<line_0_shape_set::domain_t, line_domain>::value, "Line 0 shape set domain failure");
static_assert(line_0_shape_set::num_nodes == 1, "Line 0 shape set number of nodes failure");
static_assert(line_0_shape_set::polynomial_order == 0, "Line 0 shape set polynomial order failure");
static_assert(line_0_shape_set::jacobian_order == 0, "Line 0 shape set Jacobian order failure");
static_assert(line_0_shape_set::id == 1201, "Line 0 shape set id failure");
static_assert(std::is_same<line_0_shape_set::scalar_t, double>::value, "Line 0 shape set scalar failure");
static_assert(std::is_same<line_0_shape_set::xi_t, Eigen::Matrix<double, 1, 1> >::value, "Line 0 shape set xi_t failure");
static_assert(std::is_same<line_0_shape_set::shape_t, Eigen::Matrix<double, 1, 1> >::value,
	"Line 0 shape set shape value type failure");
static_assert(std::is_same<line_0_shape_set::dshape_t, Eigen::Matrix<double, 1, 1> >::value,
	"Line 0 shape set shape derivative value type failure");
static_assert(std::is_same<line_0_shape_set::ddshape_t, Eigen::Matrix<double, 1, 1> >::value,
	"Line 0 shape set shape 2nd derivative value type failure");
static_assert(
	std::is_same<
		decltype(line_0_shape_set::template eval_shape<0>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 1, 1> const &
	>::value, "Line 0 shape set shape return type failure");
static_assert(
	std::is_same<
		decltype(line_0_shape_set::template eval_shape<1>(Eigen::Matrix<double, 1, 1>())),
		decltype(Eigen::Matrix<double, 1, 1>::Zero())
	>::value, "Line 0 shape set shape derivative return type failure");
static_assert(
	std::is_same<
		decltype(line_0_shape_set::template eval_shape<2>(Eigen::Matrix<double, 1, 1>())),
		decltype(Eigen::Matrix<double, 1, 1>::Zero())
	>::value, "Line 0 shape set shape second derivative return type failure");
static_assert(std::is_same<line_0_shape_set::position_dof_vector, tmp::vector<dof1> >::value,
		"Line 0 shape set poisition DOF vector failure");

TEST(ShapeSet, Line0){
	// test nodal location
	EXPECT_EQ(line_0_shape_set::corner_at(0)(0), 0.);
	// test shape functions in nodes
	for (unsigned i = 0; i < line_0_shape_set::num_nodes; ++i)
	{
		auto L = line_0_shape_set::eval_shape<0>(line_0_shape_set::corner_at(i));
		EXPECT_EQ(L, decltype(L)::Unit(i));
	}
	// test sum of shape functions and derivatives in random locations
	unsigned const N(1000);
	double const eps(1e-10);
	for (unsigned n = 0; n < N; ++n)
	{
		auto xi = line_0_shape_set::xi_t::Random();
		// sum of L should be one
		auto L = line_0_shape_set::eval_shape<0>(xi);
		EXPECT_NEAR(L.sum(), 1., eps);
		// sum of dL should be 0
		auto dL = line_0_shape_set::eval_shape<1>(xi);
		EXPECT_NEAR(dL.sum(), 0., eps);
		// sum of dL should be 0
		auto ddL = line_0_shape_set::eval_shape<2>(xi);
		EXPECT_NEAR(ddL.sum(), 0., eps);
	}
}



// testing line_1_shape_set
static_assert(std::is_same<line_1_shape_set::domain_t, line_domain>::value, "Line 1 shape set domain failure");
static_assert(line_1_shape_set::num_nodes == 2, "Line 1 shape set number of nodes failure");
static_assert(line_1_shape_set::polynomial_order == 1, "Line 1 shape set polynomial order failure");
static_assert(line_1_shape_set::jacobian_order == 0, "Line 1 shape set Jacobian order failure");
static_assert(line_1_shape_set::id == 1202, "Line 1 shape set id failure");
static_assert(std::is_same<line_1_shape_set::scalar_t, double>::value, "Line 1 shape set scalar failure");
static_assert(std::is_same<line_1_shape_set::xi_t, Eigen::Matrix<double, 1, 1> >::value, "Line 0 shape set xi_t failure");
static_assert(std::is_same<line_1_shape_set::shape_t, Eigen::Matrix<double, 2, 1> >::value,
	"Line 1 shape set shape value type failure");
static_assert(std::is_same<line_1_shape_set::dshape_t, Eigen::Matrix<double, 2, 1> >::value,
	"Line 1 shape set shape derivative value type failure");
static_assert(std::is_same<line_1_shape_set::ddshape_t, Eigen::Matrix<double, 2, 1> >::value,
	"Line 1 shape set shape 2nd derivative value type failure");
static_assert(
	std::is_same<
		decltype(line_1_shape_set::template eval_shape<0>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 2, 1>
	>::value, "Line 1 shape set shape return type failure");
static_assert(
	std::is_same<
		decltype(line_1_shape_set::template eval_shape<1>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 2, 1> const &
	>::value, "Line 1 shape set shape derivative return type failure");
static_assert(
	std::is_same<
		decltype(line_1_shape_set::template eval_shape<2>(Eigen::Matrix<double, 1, 1>())),
		decltype(Eigen::Matrix<double, 2, 1>::Zero())
	>::value, "Line 1 shape set shape second derivative return type failure");
static_assert(std::is_same<line_1_shape_set::position_dof_vector, tmp::vector<dof0, dof0> >::value,
		"Line 1 shape set poisition DOF vector failure");

TEST(ShapeSet, Line1){
	// test nodal location
	EXPECT_EQ(line_1_shape_set::corner_at(0)(0), -1.);
	EXPECT_EQ(line_1_shape_set::corner_at(1)(0), 1.);
	// test shape functions in nodes
	for (unsigned i = 0; i < line_1_shape_set::num_nodes; ++i)
	{
		auto L = line_1_shape_set::eval_shape<0>(line_1_shape_set::corner_at(i));
		EXPECT_EQ(L, decltype(L)::Unit(i));
	}
	// test sum of shape functions and derivatives in random locations
	unsigned const N(1000);
	double const eps(1e-10);
	for (unsigned n = 0; n < N; ++n)
	{
		auto xi = line_1_shape_set::xi_t::Random();
		auto L = line_1_shape_set::eval_shape<0>(xi);
		EXPECT_NEAR(L.sum(), 1., eps);
		auto dL = line_1_shape_set::eval_shape<1>(xi);
		EXPECT_NEAR(dL.sum(), 0., eps);
		// sum of dL should be 0
		auto ddL = line_1_shape_set::eval_shape<2>(xi);
		EXPECT_NEAR(ddL.sum(), 0., eps);
	}
}

// testing line_2_shape_set
static_assert(std::is_same<line_2_shape_set::domain_t, line_domain>::value, "Line 2 shape set domain failure");
static_assert(line_2_shape_set::num_nodes == 3, "Line 2 shape set number of nodes failure");
static_assert(line_2_shape_set::polynomial_order == 2, "Line 2 shape set polynomial order failure");
static_assert(line_2_shape_set::jacobian_order == 1, "Line 2 shape set Jacobian order failure");
static_assert(line_2_shape_set::id == 1203, "Line 2 shape set id failure");
static_assert(std::is_same<line_2_shape_set::scalar_t, double>::value, "Line 2 shape set scalar failure");
static_assert(std::is_same<line_2_shape_set::xi_t, Eigen::Matrix<double, 1, 1> >::value, "Line 0 shape set xi_t failure");
static_assert(std::is_same<line_2_shape_set::shape_t, Eigen::Matrix<double, 3, 1> >::value,
	"Line 2 shape set shape value type failure");
static_assert(std::is_same<line_2_shape_set::dshape_t, Eigen::Matrix<double, 3, 1> >::value,
	"Line 2 shape set shape derivative value type failure");
static_assert(std::is_same<line_2_shape_set::ddshape_t, Eigen::Matrix<double, 3, 1> >::value,
	"Line 2 shape set shape 2nd derivative value type failure");
static_assert(
	std::is_same<
		decltype(line_2_shape_set::template eval_shape<0>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 3, 1>
	>::value, "Line 2 shape set shape return type failure");
static_assert(
	std::is_same<
		decltype(line_2_shape_set::template eval_shape<1>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 3, 1>
	>::value, "Line 2 shape set shape derivative return type failure");
static_assert(
	std::is_same<
		decltype(line_2_shape_set::template eval_shape<2>(Eigen::Matrix<double, 1, 1>())),
		Eigen::Matrix<double, 3, 1> const &
	>::value, "Line 2 shape set shape second derivative return type failure");
static_assert(std::is_same<line_2_shape_set::position_dof_vector, tmp::vector<dof0, dof1, dof0> >::value,
		"Line 2 shape set poisition DOF vector failure");

TEST(ShapeSet, Line2){
	// test nodal location
	EXPECT_EQ(line_2_shape_set::corner_at(0)(0), -1.);
	EXPECT_EQ(line_2_shape_set::corner_at(1)(0), 0.);
	EXPECT_EQ(line_2_shape_set::corner_at(2)(0), 1.);
	// test shape functions in nodes
	for (unsigned i = 0; i < line_2_shape_set::num_nodes; ++i)
	{
		auto L = line_2_shape_set::eval_shape<0>(line_2_shape_set::corner_at(i));
		EXPECT_EQ(L, decltype(L)::Unit(i));
	}
	// test sum of shape functions and derivatives in random locations
	unsigned const N(1000);
	double const eps(1e-10);
	for (unsigned n = 0; n < N; ++n)
	{
		auto xi = line_2_shape_set::xi_t::Random();
		auto L = line_2_shape_set::eval_shape<0>(xi);
		EXPECT_NEAR(L.sum(), 1., eps);
		auto dL = line_2_shape_set::eval_shape<1>(xi);
		EXPECT_NEAR(dL.sum(), 0., eps);
		// sum of dL should be 0
		auto ddL = line_2_shape_set::eval_shape<2>(xi);
		EXPECT_NEAR(ddL.sum(), 0., eps);
	}
}

