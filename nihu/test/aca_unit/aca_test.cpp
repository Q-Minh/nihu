#include <iostream>
#include "nihu/aca/aca.hpp"

void diadic_test(void)
{
	Eigen::Matrix<double, 3, 1> a, b;
	a << 1.0, 2.0, 3.0;
	b << 1.0, 3.0, 5.0;
	Eigen::Matrix<double, 3, 3> M = a * b.transpose();
	Eigen::Matrix<double, 3, 3> U, V;

	auto R = ACA::low_rank_approx(M, 3, 3, 1e-3, 3, U, V);
	
	auto const &u = U.leftCols(R);
	auto const &v = V.leftCols(R);

	std::cout << "U:\n" << u << std::endl;
	std::cout << "V:\n" << v << std::endl;
	std::cout << "UV':\n" << u * v.transpose() << std::endl;
	std::cout << "eps: " << ((u * v.transpose()) - M).norm()/M.norm() << std::endl;
}

int main(void)
{
	diadic_test();

	return 0;
}

