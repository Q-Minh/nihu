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
			std::cout << std::endl << domain_traits::name<domain>::value << " (" << domain::id << ")\n"
                << "==========" << std::endl;
			std::cout << "domain corners: " << std::endl;
			for (unsigned i = 0; i < domain::num_corners; ++i)
				std::cout << domain::get_corner(i).transpose() << std::endl;
			std::cout << "domain center: " << domain::get_center().transpose() << std::endl;
			std::cout << "domain volume: " << domain::get_volume() << std::endl;
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
