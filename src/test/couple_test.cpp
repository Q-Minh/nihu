#include "../bem/couple.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <complex>

typedef std::complex<double> dcomplex;
typedef Eigen::Matrix<double, 3, 3> mat_t;

int main(void)
{
	mat_t m;
	dcomplex(1.0, 2.0) * m;

	return 0;
}
