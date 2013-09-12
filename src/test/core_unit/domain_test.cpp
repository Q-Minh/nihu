// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include "core/domain.hpp"

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
			std::cout << std::endl;
		}
	};
};

int main(void)
{
	std::cout << "Domains" << std::endl;
	std::cout << "=======" << std::endl;
	std::cout << "<line, tria, quad, brick>" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing domain informations" << std::endl;
	std::cout << "===========================" << std::endl;
	tmp::call_each<
		tmp::vector<line_domain, tria_domain, quad_domain, brick_domain>,
		tester<tmp::_1>
	>();
	std::cout << std::endl;
	
	return 0;
}

