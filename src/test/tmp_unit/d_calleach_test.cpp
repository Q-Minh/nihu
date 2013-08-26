#include "tmp/algorithm.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include <type_traits>
#include <iostream>

#include <iostream>

template <class T, class U>
struct tester
{
	struct type
	{
		void operator () (int const &t, int const &u)
		{
			std::cout << "T : " << T::value << "   t = " << t << 
				"     U : " << U::value << "   u = " << u << std::endl;
		}
	};
};


int main(void)
{
	typedef tmp::vector<
		std::integral_constant<int, 1>,
		std::integral_constant<int, 2>,
		std::integral_constant<int, -3>
	> vec_1;
	
	typedef tmp::vector<
		std::integral_constant<int, -4>,
		std::integral_constant<int, -6>
	> vec_2;
	
	int t = 3, u = 5;
	
	tmp::d_call_each<
		vec_1,
		vec_2,
		tester<tmp::_1, tmp::_2>,
		int,
		int
	>(t, u);
	
	return 0;
}

