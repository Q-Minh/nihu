#include "hg.hpp"

#include <Eigen/Dense>

#include <iostream>

int main(void)
{
	hg<double, int> a(1.0, 1), b(2.0, 2);

	a = a+b;

	std::cout << a <<std::endl;

	typedef Eigen::Matrix<double, 1, 2> vec_t;
	hg<vec_t, vec_t> v;

	return 0;
}
