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

#include "nihu/util/conditional_precompute.hpp"

#include <iostream>
#include <Eigen/Core>

typedef Eigen::Matrix<double, 3, 3> dMatrix;

struct Func
{
	static dMatrix eval(char arg)
	{
		return arg == 'a' ? dMatrix::Zero() : dMatrix::Ones();
	}
};

typedef NiHu::conditional_precompute<NiHu::matrix_function_complexity::general, Func, char> fly_t;
typedef NiHu::conditional_precompute<NiHu::matrix_function_complexity::constant, Func, char> store_t;

typedef NiHu::conditional_precompute_instance<NiHu::matrix_function_complexity::general, Func, char> fly_inst_t;
typedef NiHu::conditional_precompute_instance<NiHu::matrix_function_complexity::constant, Func, char> store_inst_t;

int main(void)
{
	std::cout << fly_t::eval('b') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	std::cout << fly_t::eval('a') << std::endl;
	std::cout << store_t::eval('b') << std::endl;

	fly_inst_t fly_inst;
	std::cout << fly_inst('a') << std::endl;
	std::cout << fly_inst('b') << std::endl;

	store_inst_t store_inst('a');
	std::cout << store_inst() << std::endl;

	return 0;
}
