#include "../bem/couple.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <complex>

typedef std::complex<double> dcomplex;
typedef Eigen::Matrix<double, 3, 3> mat_t;

int main(void)
{
	/*
	couple<dcomplex> a(dcomplex(1.0, 2.0), dcomplex(3.0, 4.0));
	mat_t m = mat_t::Ones();
	std::cout << (m * a * m).first();
	*/
	return 0;
}
