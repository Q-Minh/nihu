#include "core/space.hpp"
#include <gtest/gtest.h>

// testing if space dimensions are defined correctly
static_assert(space<double, 1>::dimension == 1, "Space1D dimension error");
static_assert(space<double, 2>::dimension == 2, "Space2D dimension error");
static_assert(space<double, 3>::dimension == 3, "Space3D dimension error");
// testing if space scalar is defined correctly
static_assert(std::is_same<space<double, 1>::scalar_t, double>::value, "Space scalar error");
static_assert(std::is_same<space<float, 1>::scalar_t, float>::value, "Space scalar error");
static_assert(std::is_same<space<int, 1>::scalar_t, int>::value, "Space scalar error");
// testing if space shorthands are defined correctly
static_assert(std::is_same<space_1d<float>, space<float, 1> >::value, "Space_1d shorthand error");
static_assert(std::is_same<space_2d<float>, space<float, 2> >::value, "Space_2d shorthand error");
static_assert(std::is_same<space_3d<float>, space<float, 3> >::value, "Space_3d shorthand error");
// testing if default space shorthands are defined correctly
static_assert(std::is_same<space_1d<>, space<double, 1> >::value, "Space_1d default shorthand error");
static_assert(std::is_same<space_2d<>, space<double, 2> >::value, "Space_2d default shorthand error");
static_assert(std::is_same<space_3d<>, space<double, 3> >::value, "Space_3d default shorthand error");

#include "library/lib_domain.hpp"

// testing line_1_domain
static_assert(std::is_same<line_domain::space_t, space_1d<> >::value, "Line domain space failure");
static_assert(line_domain::num_corners == 2, "Line domain number of corners failure");
static_assert(line_domain::id == 12, "Line domain ID failure");
static_assert(line_domain::get_id() == 12, "Line domain get_ID() failure");
static_assert(line_domain::dimension == 1, "Line domain dimension failure");
static_assert(std::is_same<line_domain::scalar_t, double>::value, "Line domain scalar failure");
static_assert(std::is_same<line_domain::xi_t, space_1d<>::location_t>::value, "Line domain location failure");
static_assert(line_domain::get_volume()==2., "Line domain get_volume() failure");

TEST(Domain, Line){
	EXPECT_EQ(line_domain::get_center()(0), 0.);
	EXPECT_EQ(line_domain::get_corners()[0](0), -1.);
	EXPECT_EQ(line_domain::get_corners()[1](0), +1.);
	EXPECT_EQ(line_domain::get_corner(0)(0), -1.);
	EXPECT_EQ(line_domain::get_corner(1)(0), +1.);
}

// testing tria_domain
static_assert(std::is_same<tria_domain::space_t, space_2d<> >::value, "Tria domain space failure");
static_assert(tria_domain::num_corners == 3, "Tria domain number of corners failure");
static_assert(tria_domain::id == 23, "Tria domain ID failure");
static_assert(tria_domain::get_id() == 23, "Tria domain get_ID() failure");
static_assert(tria_domain::dimension == 2, "Tria domain dimension failure");
static_assert(std::is_same<tria_domain::scalar_t, double>::value, "Tria domain scalar failure");
static_assert(std::is_same<tria_domain::xi_t, space_2d<>::location_t>::value, "Tria domain location failure");
static_assert(tria_domain::get_volume()==.5, "Tria domain get_volume() failure");

TEST(Domain, Tria){
	EXPECT_EQ(tria_domain::get_center()(0), 1./3.);
	EXPECT_EQ(tria_domain::get_center()(1), 1./3.);
	tria_domain::scalar_t corners[][2] = { {0., 0.}, { 1., 0.}, { 0.,  1.} };
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 2; ++j)
			EXPECT_EQ(tria_domain::get_corner(i)(j), corners[i][j]);
}

// testing quad_domain
static_assert(std::is_same<quad_domain::space_t, space_2d<> >::value, "Quad domain space failure");
static_assert(quad_domain::num_corners == 4, "Quad domain number of corners failure");
static_assert(quad_domain::id == 24, "Quad domain ID failure");
static_assert(quad_domain::get_id() == 24, "Quad domain get_ID() failure");
static_assert(quad_domain::dimension == 2, "Quad domain dimension failure");
static_assert(std::is_same<quad_domain::scalar_t, double>::value, "Quad domain scalar failure");
static_assert(std::is_same<quad_domain::xi_t, space_2d<>::location_t>::value, "Quad domain location failure");
static_assert(quad_domain::get_volume()==4., "Quad domain get_volume() failure");

TEST(Domain, Quad){
	EXPECT_EQ(quad_domain::get_center()(0), 0.);
	EXPECT_EQ(quad_domain::get_center()(1), 0.);
	quad_domain::scalar_t corners[][2] = { {-1., -1.}, { 1., -1.}, { 1.,  1.}, {-1.,  1.} };
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 2; ++j)
			EXPECT_EQ(quad_domain::get_corner(i)(j), corners[i][j]);
}

// testing brick_domain
static_assert(std::is_same<brick_domain::space_t, space_3d<> >::value, "Brick domain space failure");
static_assert(brick_domain::num_corners == 8, "Brick domain number of corners failure");
static_assert(brick_domain::id == 38, "Brick domain ID failure");
static_assert(brick_domain::get_id() == 38, "Brick domain get_ID() failure");
static_assert(brick_domain::dimension == 3, "Brick domain dimension failure");
static_assert(std::is_same<brick_domain::scalar_t, double>::value, "Brick domain scalar failure");
static_assert(std::is_same<brick_domain::xi_t, space_3d<>::location_t>::value, "Brick domain location failure");
static_assert(brick_domain::get_volume()==8., "Brick domain get_volume() failure");

TEST(Domain, Brick){
	EXPECT_EQ(brick_domain::get_center()(0), 0.);
	EXPECT_EQ(brick_domain::get_center()(1), 0.);
	brick_domain::scalar_t corners[][3] = {
	{-1., -1., -1.},
	{ 1., -1., -1.},
	{ 1.,  1., -1.},
	{-1.,  1., -1.},
	{-1., -1.,  1.},
	{ 1., -1.,  1.},
	{ 1.,  1.,  1.},
	{-1.,  1.,  1.}
	};
	for (unsigned i = 0; i < 8; ++i)
		for (unsigned j = 0; j < 3; ++j)
			EXPECT_EQ(brick_domain::get_corner(i)(j), corners[i][j]);
}

