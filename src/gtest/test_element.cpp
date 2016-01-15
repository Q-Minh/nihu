/// \file test_element.cpp
/// Gtest file of class Element

#include <gtest/gtest.h>
#include "library/lib_element.hpp"

template <class ElemType>
struct Stored;

template <>
struct Stored<line_1_elem>
{
	Stored(void)
	{
		coords << 1.0, 2.0, 2.0, 4.0;
		center << 1.5, 3.0;
		size_estimate = std::sqrt(5.0);
		xi << .3;
		x << ((1.0-xi(0)) * 1.0 + (1.0+xi(0)) * 2.0) / 2.0, ((1.0-xi(0)) * 2.0 + (1.0+xi(0)) * 4.0) / 2.0;
		normal << 2.0/2.0, -1.0/2.0;
	}

	line_1_elem::coords_t coords;
	line_1_elem::x_t center;
	line_1_elem::xi_t xi;
	line_1_elem::x_t x;
	line_1_elem::x_t normal;
	double size_estimate;
};


template <class ElemType>
void tester(void)
{
	Stored<ElemType> stored;

	ElemType e(stored.coords);

	double eps = 1e-8;

	EXPECT_NEAR((e.get_center() - stored.center).norm(), 0.0, eps);
	EXPECT_NEAR(e.get_linear_size_estimate(), stored.size_estimate, eps);
	EXPECT_NEAR((e.get_x(stored.xi) - stored.x).norm(), 0.0, eps);
}

TEST(Element, Line1Elem)
{
	tester<line_1_elem>();
}

