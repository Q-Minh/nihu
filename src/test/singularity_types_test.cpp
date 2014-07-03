#include <iostream>
#include <core/singularity_types.hpp>

int main(void)
{
	std::cout << tmp::greater<
		singularity_type::log<1>,
		singularity_type::log<2>
	>::value << std::endl;

	std::cout << tmp::less<
		singularity_type::log<1>,
		singularity_type::log<2>
	>::value << std::endl;

	typedef tmp::max_<
		singularity_type::log<1>, singularity_type::log<2>, singularity_type::inverse<3>
	>::type max_sing_t;

	std::cout << std::is_same<
		max_sing_t, singularity_type::inverse<3>
	>::value << std::endl;

	typedef tmp::min_<
		singularity_type::log<1>, singularity_type::log<2>, singularity_type::inverse<3>
	>::type min_sing_t;

	std::cout << std::is_same<
		min_sing_t, singularity_type::log<1>
	>::value << std::endl;

	return 0;
}
