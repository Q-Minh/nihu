#include "util/conditional_precompute.hpp"

#include <iostream>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, 3, 3> dMatrix;

struct Func
{
	static dMatrix eval(char arg)
	{
		return arg == 'a' ? dMatrix::Zero() : dMatrix::Ones();
	}
};

typedef conditional_precompute<matrix_function_complexity::general, Func, char> fly_t;
typedef conditional_precompute<matrix_function_complexity::constant, Func, char> store_t;

typedef conditional_precompute_instance<matrix_function_complexity::general, Func, char> fly_inst_t;
typedef conditional_precompute_instance<matrix_function_complexity::constant, Func, char> store_inst_t;

int main(void)
{
	std::cout << fly_t::eval('b') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	std::cout << fly_t::eval('a') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	fly_inst_t fly_inst;
	std::cout << fly_inst('a') << std::endl;
	std::cout << fly_inst('b') << std::endl;

	store_inst_t store_inst('a');
	std::cout << store_inst() << std::endl;

	return 0;
}
