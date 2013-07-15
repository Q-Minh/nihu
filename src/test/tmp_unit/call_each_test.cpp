#include "tmp/algorithm.hpp"
#include "tmp/control.hpp"
#include <type_traits>
#include <iostream>

template <class T>
struct tester
{
	struct type
	{
		void operator () (int const &a, char const &b)
		{
			std::cout << T::value << a << b << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<
		std::integral_constant<int, 1>,
		std::integral_constant<int, 2>
	> vec;

	typedef tester<tmp::_1> test_t;
	
	tmp::call_each<
		vec,
		test_t,
		int const &,
		char const &
	>(1, 'a');
}
