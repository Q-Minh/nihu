#include <gtest/gtest.h>

#include "nihu/fmm/bounding_box.hpp"

template <size_t Dim>
void bb_tester()
{
	static size_t const dimension = Dim;
	typedef NiHu::fmm::bounding_box<dimension, double> bounding_box_t;

	// default constructor
	bounding_box_t bb0;

	// 2param constructor
	bounding_box_t bb1(bounding_box_t::location_t::Constant(0.0), 2.0);

	// check if these boxes are adjacent
	EXPECT_TRUE(bb0.is_adjacent(bb1));
	EXPECT_TRUE(bb1.is_adjacent(bb0));

	// constructor from nodes
	size_t n_nodes = 10;
	Eigen::Matrix<double, dimension, Eigen::Dynamic> nodes(dimension, n_nodes);
	nodes.setZero();
	nodes.row(0).setLinSpaced(n_nodes, -1., 1.);
	bounding_box_t bb3(nodes);
	EXPECT_DOUBLE_EQ(2.0, bb3.get_diameter());

	// copy constructor
	bounding_box_t bb2 = bb1;

	// assignment operator
	bb2 = bb1;
}

TEST(bounding_box, _1d)
{
	bb_tester<1>();
}

TEST(bounding_box, _2d)
{
	bb_tester<2>();
}

TEST(bounding_box, _3d)
{
	bb_tester<3>();
}
