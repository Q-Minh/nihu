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

#include "nihu/tmp/interval.hpp"
#include "nihu/tmp/control.hpp"
#include <iostream>

typedef tmp::vector<
	tmp::break_point<std::ratio<2>, tmp::integer<int, 8> >,
	tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 1> >
> inter1_t;

typedef tmp::vector<
	tmp::break_point<std::ratio<3>, tmp::integer<int, 6> >,
	tmp::break_point<std::ratio<5>, tmp::integer<int, 4> >,
	tmp::break_point<tmp::ratio_infinite, tmp::integer<int, 2> >
> inter2_t;

typedef tmp::merge_intervals<inter1_t, inter2_t>::type merged_t;

template <class BP>
struct print_break_point
{
	struct type
	{
		void operator()(void)
		{
			std::cout << BP::x_value() << " : " << BP::y::value << std::endl;
		}
	};
};

int main(void)
{
	std::cout << "Test intervals" << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "I1 :" << std::endl;
	tmp::call_each<inter1_t, print_break_point<tmp::_1> >();
	std::cout << "I2 :" << std::endl;
	tmp::call_each<inter2_t, print_break_point<tmp::_1> >();
	std::cout << std::endl;

	std::cout << "Testing interval merge" << std::endl;
	std::cout << "======================" << std::endl;
	
	std::cout << "Merged : " << std::endl;
	tmp::call_each<merged_t, print_break_point<tmp::_1> >();
	std::cout << std::endl;

	std::cout << "Testing interval evaluation" << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << "  x : I1, I2, Merged " << std::endl;
	std::cout <<  .2 << " : " <<' ' << tmp::eval_interval<inter1_t>( .2) << ' ' << tmp::eval_interval<inter2_t>( .2) << ' ' << tmp::eval_interval<merged_t>( .2) << std::endl;
	std::cout << 2.2 << " : " <<' ' << tmp::eval_interval<inter1_t>(2.2) << ' ' << tmp::eval_interval<inter2_t>(2.2) << ' ' << tmp::eval_interval<merged_t>(2.2) << std::endl;
	std::cout << 4.2 << " : " <<' ' << tmp::eval_interval<inter1_t>(4.2) << ' ' << tmp::eval_interval<inter2_t>(4.2) << ' ' << tmp::eval_interval<merged_t>(4.2) << std::endl;
	std::cout << 9.8 << " : " <<' ' << tmp::eval_interval<inter1_t>(9.8) << ' ' << tmp::eval_interval<inter2_t>(9.8) << ' ' << tmp::eval_interval<merged_t>(9.8) << std::endl;
	std::cout << std::endl;

	return 0;
}
