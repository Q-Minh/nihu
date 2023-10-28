/// \file test_element.cpp
/// Gtest file of class Element

#include <gtest/gtest.h>
#include "library/lib_element.hpp"

template <class ElemType>
struct Stored;

template <>
struct Stored<NiHu::line_1_elem>
{
	Stored(void)
	{
		coords <<
			1.0, 2.0,
			2.0, 4.0;
		center << 1.5, 3.0;
		size_estimate = std::sqrt(5.0);
		xi << .3;
		x << ((1.0-xi(0)) * 1.0 + (1.0+xi(0)) * 2.0) / 2.0, ((1.0-xi(0)) * 2.0 + (1.0+xi(0)) * 4.0) / 2.0;
		normal << 2.0/2.0, -1.0/2.0;
	}

	NiHu::line_1_elem::coords_t coords;
	NiHu::line_1_elem::x_t center;
	NiHu::line_1_elem::xi_t xi;
	NiHu::line_1_elem::x_t x;
	NiHu::line_1_elem::x_t normal;
	double size_estimate;
};


template <>
struct Stored<NiHu::line_2_elem>
{
	Stored(void)
	{
		coords <<
			0.0, 1.0, 2.0,
			0.0, 1.0, 1.0;
		center << 1.0, 1.0;
		size_estimate = std::sqrt(5.0);
		xi << .3;
		x << 1.3, 1.105;
		normal << .2, -1.0;
	}

	NiHu::line_2_elem::coords_t coords;
	NiHu::line_2_elem::x_t center;
	NiHu::line_2_elem::xi_t xi;
	NiHu::line_2_elem::x_t x;
	NiHu::line_2_elem::x_t normal;
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
	tester<NiHu::line_1_elem>();
}

TEST(Element, Line2Elem)
{
	tester<NiHu::line_2_elem>();
}

