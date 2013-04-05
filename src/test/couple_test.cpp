#include "../bem/couple.hpp"

#include <iostream>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, 3, 3> mat_t;

int main(void)
{
	couple<mat_t> a(mat_t::Zero(), mat_t::Ones()), b(mat_t::Zero(), mat_t::Ones());
	std::cout << (a*mat_t::Ones()).second() << std::endl;
	std::cout << (mat_t::Ones()*a).second() << std::endl;
	return 0;
}
