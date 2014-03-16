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

#include <iostream>
#include "util/couple.hpp"
#include <Eigen/Dense>

int main(void)
{
	couple<int, float, char, Eigen::Matrix<double, 3, 3> > a;
	std::cout << a.get<0>() << '\n';
	std::cout << a.get<1>() << '\n';
	std::cout << a.get<2>() << '\n';
	std::cout << a.get<3>() << '\n';

	auto b = a * 2;

	std::cout << b.get<0>() << '\n';
	std::cout << b.get<1>() << '\n';
	std::cout << b.get<2>() << '\n';
	std::cout << b.get<3>() << '\n';
}

