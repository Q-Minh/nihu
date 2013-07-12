#include "tmp/control.hpp"
#include "bem/space.hpp"

#include <iostream>

template <class space>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << space::dimension << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<space_2d, space_3d>,
		tester<tmp::_1>
	>();
	
	return 0;
}

