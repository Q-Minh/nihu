/// \file test_cluster.cpp
/// \brief Gtest file for class NiHu::fmm::cluster

#include <gtest/gtest.h>

#include "fmm/empty_cluster.hpp"

TEST(cluster, empty)
{
	NiHu::fmm::empty_cluster<1> c[3];
	c[0].set_level(0);
	c[1].set_level(1);
	c[2].set_level(1);

	// check parent-child relations
	c[1].set_parent(0);
	c[2].set_parent(0);
	c[0].add_child(1);
	c[0].add_child(2);
	// cannot add a child twice
	EXPECT_THROW(c[0].add_child(1), std::logic_error);

	// check roots and leaves
	EXPECT_TRUE(c[0].is_root());
	EXPECT_FALSE(c[1].is_root());
	EXPECT_FALSE(c[2].is_root());
	EXPECT_FALSE(c[0].is_leaf());
	EXPECT_TRUE(c[1].is_leaf());
	EXPECT_TRUE(c[2].is_leaf());

	// check source and receiver clusters
	// initially no cluster is source or receiver
	for (size_t i = 0; i < 3; ++i)
	{
		EXPECT_FALSE(c[i].is_source());
		EXPECT_FALSE(c[i].is_receiver());
	}
	std::vector<size_t> idx(20);
	for (size_t i = 0; i < 20; ++i)
		idx.push_back(i);
	// c[0] contains everybody
	c[0].set_src_node_idx(idx);
	c[0].set_rec_node_idx(idx);
	// c[1] contains first half
	c[1].set_src_node_idx(idx.begin(), idx.begin() + 10);
	c[1].set_rec_node_idx(idx.begin(), idx.begin() + 10);
	// c[2] contains second half
	c[2].set_src_node_idx(idx.begin()+10, idx.end());
	c[2].set_rec_node_idx(idx.begin()+10, idx.end());
	for (size_t i = 0; i < 3; ++i)
	{
		EXPECT_TRUE(c[i].is_source());
		EXPECT_TRUE(c[i].is_receiver());
	}
}
