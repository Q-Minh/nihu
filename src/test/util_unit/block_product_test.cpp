#include "util/block_product.hpp"
#include <iostream>

int main(void)
{
	Eigen::Matrix<double, 2, 1> v1;
	Eigen::Matrix<double, 3, 3> m;
	Eigen::Matrix<double, 4, 1> v2;

	std::cout << block_product(v1, m, v2) << std::endl;

	Eigen::Matrix<double, 1, 1> vv1;
	vv1 << 1.0;
	double mat = 2.0;
	Eigen::Matrix<double, 1, 1> vv2;
	vv2 << 1.0;

	std::cout << block_product(vv1, mat, vv2) << std::endl;
}
