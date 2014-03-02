#include "util/block_product.hpp"
#include <iostream>

int main(void)
{
	Eigen::Matrix<double, 2, 1> v1;
	Eigen::Matrix<double, 3, 3> m;
	Eigen::Matrix<double, 4, 1> v2;

	std::cout << block_product(v1, m, v2) << std::endl;
}
