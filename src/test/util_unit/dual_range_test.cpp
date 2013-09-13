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

#include "util/dual_range.hpp"
#include <iostream>

template <class range>
void iterate(range r)
{
	for (auto it = r.begin(); it != r.end(); ++it)
		std::cout << *(it.get_first()) << ' ' << *(it.get_second()) << std::endl;
}

int main(void)
{
	int a[] = {0, 1, 2, 3};
	char b[] = {'a', 'b', 'c', 'd'};
	
	std::cout << "diagonal iteration" << std::endl;
	iterate(create_dual_range(iteration::diagonal(), a, a+4, b, b+4));

	std::cout << "matrix iteration" << std::endl;
	iterate(create_dual_range(iteration::diadic(), a, a+4, b, b+4));

	std::cout << "matrix iteration with empty inner" << std::endl;
	iterate(create_dual_range(iteration::diadic(), a, a+4, b, b));

	std::cout << "matrix iteration with empty outer" << std::endl;
	iterate(create_dual_range(iteration::diadic(), a, a, b, b+4));
	
	return 0;
}

