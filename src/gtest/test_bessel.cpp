#include "util/math_functions.hpp"


#include <gtest/gtest.h>

void bessel_tester(double data[][4], unsigned size, std::complex<double> (*func)(std::complex<double> const &) )
{
	for (unsigned i = 0; i < size; ++i)
	{
		double const eps(1e-6);
		std::complex<double> z(data[i][0], data[i][1]);
		std::complex<double> Mat(data[i][2], data[i][3]);
		std::complex<double> NiHu = func(z);
		double E = std::abs((Mat-NiHu)/Mat);
		if (E > eps)
			std::cout << "R: " << std::abs(z) << " phi: " << std::arg(z)/M_PI*180. << " z: " << z << " Mat: " << Mat << " NiHu: " << NiHu << " eps " << E << '\n';
		EXPECT_NEAR(E, 0., eps);
	}
}

TEST(Bessel, J0) {
	#include "data/besselJ0.dat"
	bessel_tester(besselJ0, sizeof(besselJ0)/sizeof(besselJ0[0]), bessel::J<0>);
}

TEST(Bessel, J1) {
	#include "data/besselJ1.dat"
	bessel_tester(besselJ1, sizeof(besselJ1)/sizeof(besselJ1[0]), bessel::J<1>);
}

TEST(Bessel, Y0) {
	#include "data/besselY0.dat"
	bessel_tester(besselY0, sizeof(besselY0)/sizeof(besselY0[0]), bessel::Y<0>);
}

TEST(Bessel, Y1) {
	#include "data/besselY1.dat"
	bessel_tester(besselY1, sizeof(besselY1)/sizeof(besselY1[0]), bessel::Y<1>);
}

TEST(Bessel, H01) {
	#include "data/besselH01.dat"
	bessel_tester(besselH01, sizeof(besselH01)/sizeof(besselH01[0]), bessel::H<0,1>);
}

TEST(Bessel, H02) {
	#include "data/besselH02.dat"
	bessel_tester(besselH02, sizeof(besselH02)/sizeof(besselH02[0]), bessel::H<0,2>);
}

TEST(Bessel, H11) {
	#include "data/besselH11.dat"
	bessel_tester(besselH11, sizeof(besselH11)/sizeof(besselH11[0]), bessel::H<1,1>);
}

TEST(Bessel, H12) {
	#include "data/besselH12.dat"
	bessel_tester(besselH12, sizeof(besselH12)/sizeof(besselH12[0]), bessel::H<1,2>);
}

TEST(Bessel, K0) {
	#include "data/besselK0.dat"
	bessel_tester(besselK0, sizeof(besselK0)/sizeof(besselK0[0]), bessel::K<0>);
}

TEST(Bessel, K1) {
	#include "data/besselK1.dat"
	bessel_tester(besselK1, sizeof(besselK1)/sizeof(besselK1[0]), bessel::K<1>);
}
