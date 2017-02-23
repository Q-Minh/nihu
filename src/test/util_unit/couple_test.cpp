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
	// testing couple constructor and get functions
	NiHu::couple<int, float, char, Eigen::Matrix<double, 3, 3> > a(3, 3.14, 'b', Eigen::Matrix<double, 3, 3>::Constant(2.0));
	std::cout << a << std::endl;
	
	// testing eval_to_tuple
	auto t = a.eval_to_tuple();
	std::cout << std::get<0>(t) << std::endl;
	std::cout << std::get<1>(t) << std::endl;
	std::cout << std::get<2>(t) << std::endl;
	std::cout << std::get<3>(t) << std::endl;

	// testing right product
	auto b = a * 2;
	std::cout << b << std::endl;
	
	// testing left product
	auto c = 2 * a;
	std::cout << c << std::endl;

	// testing couple block
	Eigen::Matrix<double, 3, 3> m1 = Eigen::Matrix<double, 3, 3>::Zero(), m2 = Eigen::Matrix<double, 3, 3>::Zero();
	auto C = NiHu::create_couple(m1, m2);
	C.block<2,2>(1,1) += NiHu::create_couple(Eigen::Matrix<double,2,2>::Constant(1.0), Eigen::Matrix<double, 2, 2>::Constant(1.0));
	std::cout << C << std::endl;
}

