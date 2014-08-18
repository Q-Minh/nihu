// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
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

#include "space_test.h"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include "core/space.hpp"

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

int space_test(void)
{
	std::cout << "Spaces" << std::endl;
	std::cout << "======" << std::endl;
	std::cout << "<space_1d, space_2d, space_3d>" << std::endl;
	std::cout << std::endl;

	std::cout << "Testing space dims and loc types" << std::endl;
	std::cout << "================================" << std::endl;

	tmp::call_each<
		tmp::vector<space_1d<>, space_2d<>, space_3d<> >,
		tester<tmp::_1>
	>();
	std::cout << std::endl;

	return 0;
}

