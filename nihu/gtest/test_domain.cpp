/// \file test_domain
/// Gtest file of class NiHu::domain

#include <gtest/gtest.h>

#include "nihu/library/lib_domain.hpp"

// testing line_domain
static_assert(std::is_same<NiHu::line_domain::space_t,
	NiHu::space_1d<> >::value, "Line domain space failure");
static_assert(NiHu::line_domain::num_corners == 2,
	"Line domain number of corners failure");
static_assert(NiHu::line_domain::id == 12, "Line domain ID failure");
static_assert(NiHu::line_domain::get_id() == 12,
	"Line domain get_ID() failure");
static_assert(NiHu::line_domain::dimension == 1,
	"Line domain dimension failure");
static_assert(std::is_same<NiHu::line_domain::scalar_t, double>::value,
	"Line domain scalar failure");
static_assert(std::is_same<NiHu::line_domain::xi_t,
	NiHu::space_1d<>::location_t>::value, "Line domain location failure");
static_assert(NiHu::line_domain::get_volume()==2.,
	"Line domain get_volume() failure");

TEST(Domain, Line){
	EXPECT_EQ(NiHu::line_domain::get_center()(0), 0.);
	EXPECT_EQ(NiHu::line_domain::get_corners()[0](0), -1.);
	EXPECT_EQ(NiHu::line_domain::get_corners()[1](0), +1.);
	EXPECT_EQ(NiHu::line_domain::get_corner(0)(0), -1.);
	EXPECT_EQ(NiHu::line_domain::get_corner(1)(0), +1.);
}

// testing tria_domain
static_assert(std::is_same<NiHu::tria_domain::space_t,
	NiHu::space_2d<> >::value, "Tria domain space failure");
static_assert(NiHu::tria_domain::num_corners == 3,
	"Tria domain number of corners failure");
static_assert(NiHu::tria_domain::id == 23, "Tria domain ID failure");
static_assert(NiHu::tria_domain::get_id() == 23,
	"Tria domain get_ID() failure");
static_assert(NiHu::tria_domain::dimension == 2,
	"Tria domain dimension failure");
static_assert(std::is_same<NiHu::tria_domain::scalar_t, double>::value,
	"Tria domain scalar failure");
static_assert(std::is_same<NiHu::tria_domain::xi_t,
	NiHu::space_2d<>::location_t>::value, "Tria domain location failure");
static_assert(NiHu::tria_domain::get_volume()==.5,
	"Tria domain get_volume() failure");

TEST(Domain, Tria){
	EXPECT_EQ(NiHu::tria_domain::get_center()(0), 1./3.);
	EXPECT_EQ(NiHu::tria_domain::get_center()(1), 1./3.);
	NiHu::tria_domain::scalar_t corners[][2] =	{
		{0., 0.}, { 1., 0.}, { 0.,  1.}
	};
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 2; ++j)
			EXPECT_EQ(NiHu::tria_domain::get_corner(i)(j), corners[i][j]);
}

// testing quad_domain
static_assert(std::is_same<NiHu::quad_domain::space_t,
	NiHu::space_2d<> >::value, "Quad domain space failure");
static_assert(NiHu::quad_domain::num_corners == 4,
	"Quad domain number of corners failure");
static_assert(NiHu::quad_domain::id == 24, "Quad domain ID failure");
static_assert(NiHu::quad_domain::get_id() == 24,
	"Quad domain get_ID() failure");
static_assert(NiHu::quad_domain::dimension == 2,
	"Quad domain dimension failure");
static_assert(std::is_same<NiHu::quad_domain::scalar_t, double>::value,
	"Quad domain scalar failure");
static_assert(std::is_same<NiHu::quad_domain::xi_t,
	NiHu::space_2d<>::location_t>::value, "Quad domain location failure");
static_assert(NiHu::quad_domain::get_volume()==4.,
	"Quad domain get_volume() failure");

TEST(Domain, Quad){
	EXPECT_EQ(NiHu::quad_domain::get_center()(0), 0.);
	EXPECT_EQ(NiHu::quad_domain::get_center()(1), 0.);
	NiHu::quad_domain::scalar_t corners[][2] = {
		{-1., -1.}, { 1., -1.}, { 1.,  1.}, {-1.,  1.}
	};
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 2; ++j)
			EXPECT_EQ(NiHu::quad_domain::get_corner(i)(j), corners[i][j]);
}

// testing brick_domain
static_assert(std::is_same<NiHu::brick_domain::space_t,
	NiHu::space_3d<> >::value, "Brick domain space failure");
static_assert(NiHu::brick_domain::num_corners == 8,
	"Brick domain number of corners failure");
static_assert(NiHu::brick_domain::id == 38, "Brick domain ID failure");
static_assert(NiHu::brick_domain::get_id() == 38,
	"Brick domain get_ID() failure");
static_assert(NiHu::brick_domain::dimension == 3, 
	"Brick domain dimension failure");
static_assert(std::is_same<NiHu::brick_domain::scalar_t, double>::value,
	"Brick domain scalar failure");
static_assert(std::is_same<NiHu::brick_domain::xi_t,
	NiHu::space_3d<>::location_t>::value, "Brick domain location failure");
static_assert(NiHu::brick_domain::get_volume()==8.,
	"Brick domain get_volume() failure");

TEST(Domain, Brick){
	EXPECT_EQ(NiHu::brick_domain::get_center()(0), 0.);
	EXPECT_EQ(NiHu::brick_domain::get_center()(1), 0.);
	NiHu::brick_domain::scalar_t corners[][3] = {
		{-1., -1., -1.}, { 1., -1., -1.}, { 1.,  1., -1.}, {-1.,  1., -1.},
		{-1., -1.,  1.}, { 1., -1.,  1.}, { 1.,  1.,  1.}, {-1.,  1.,  1.}
	};
	for (unsigned i = 0; i < 8; ++i)
		for (unsigned j = 0; j < 3; ++j)
			EXPECT_EQ(NiHu::brick_domain::get_corner(i)(j), corners[i][j]);
}

