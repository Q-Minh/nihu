#include <iostream>
#include "util/couple.hpp"
#include <Eigen/Dense>

int main(void)
{
	couple<int, float, char, Eigen::Matrix<double, 3, 3> > a;
	std::cout << a.get<0>() << '\n';
	std::cout << a.get<1>() << '\n';
	std::cout << a.get<2>() << '\n';
	std::cout << a.get<3>() << '\n';

	auto b = a * 2;

	std::cout << b.get<0>() << '\n';
	std::cout << b.get<1>() << '\n';
	std::cout << b.get<2>() << '\n';
	std::cout << b.get<3>() << '\n';
}

