#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include "bem/space.hpp"

#include <iostream>

template <class space>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << "Dim: " << space::dimension << ", ";
			std::cout << "Loc: " << space::location_t::RowsAtCompileTime << " x " << space::location_t::ColsAtCompileTime << std::endl;
		}
	};
};

int main(void)
{
	std::cout << "Spaces" << std::endl;
	std::cout << "======" << std::endl;
	std::cout << "<space_1d, space_2d, space_3d>" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing space dims and loc types" << std::endl;
	std::cout << "================================" << std::endl;

	tmp::call_each<
		tmp::vector<space_1d, space_2d, space_3d>,
		tester<tmp::_1>
	>();
	std::cout << std::endl;
	
	return 0;
}

