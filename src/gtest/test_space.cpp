/// \file test_space.cpp
/// GTest tester for class space

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


