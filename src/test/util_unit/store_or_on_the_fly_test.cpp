#include "util/conditional_precompute.hpp"

#include <iostream>

struct Func
{
	static bool eval(char arg)
	{
		return arg == 'a';
	}
};

typedef conditional_precompute<true, Func, char> fly_t;
typedef conditional_precompute<false, Func, char> store_t;

int main(void)
{
	std::cout << fly_t::eval('b') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	std::cout << fly_t::eval('a') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	return 0;
}
