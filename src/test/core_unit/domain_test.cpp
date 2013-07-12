#include "tmp/control.hpp"
#include "bem/domain.hpp"

#include <iostream>


template <class domain>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl << "domain id: " << domain::id << std::endl;
			std::cout << "domain corners: " << std::endl;
			for (unsigned i = 0; i < domain::num_corners; ++i)
				std::cout << domain::get_corner(i).transpose() << std::endl;
			std::cout << "domain center: " << std::endl;
			std::cout << domain::get_center().transpose() << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<line_domain, tria_domain, quad_domain, brick_domain>,
		tester<tmp::_1>
	>();
	
	return 0;
}

