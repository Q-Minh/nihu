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

#include "tmp/algorithm.hpp"
#include "tmp/integer.hpp"
#include "tmp/vector.hpp"
#include "tmp/print.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::vector<
		tmp::integer<int, 1>,
		tmp::integer<int, 3>,
		tmp::integer<int, 28>,
		tmp::integer<int, -32>,
		tmp::integer<int, 3>,
		tmp::integer<int, 0>,
		tmp::integer<int, -2>
	> v1;
	
	typedef tmp::vector<
		tmp::integer<int, -8>,
		tmp::integer<int, 6>,
		tmp::integer<int, 2>,
		tmp::integer<int, -1>
	> v2;
	
	
	std::cout << "Testing accumulation" << std::endl;
	std::cout << "====================" << std::endl;
	std::cout << "v1       = "; print<v1>::eval();
	std::cout << "v2       = "; print<v2>::eval();
	std::cout << "sum(v1)  = "; std::cout << tmp::accumulate<v1, tmp::integer<int, 0>>::type::value << std::endl;
	std::cout << "prod(v2) = "; std::cout << tmp::accumulate<v2, tmp::integer<int, 1>, tmp::mul<tmp::_1, tmp::_2> >::type::value << std::endl;
	std::cout << "max(v1)  = "; std::cout << tmp::max<v1>::type::value << std::endl;
	std::cout << "min(v1)  = "; std::cout << tmp::min<v1>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing bubble_sort" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "v1   = "; print<v1>::eval();
	std::cout << "v1 < = "; print<tmp::bubble_sort<v1>::type>::eval();
	std::cout << "v1 > = ";print<tmp::bubble_sort<v1, tmp::less<tmp::_1, tmp::_2> >::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing concatenate" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "v1       = "; print<v1>::eval();
	std::cout << "v2       = "; print<v2>::eval();
	std::cout << "[v2, v1] = "; print<tmp::concatenate<v2, v1>::type>::eval();
	
	std::cout << "Testing find" << std::endl;
	std::cout << "============" << std::endl; 
	std::cout << "v1            = "; print<v1>::eval();
	std::cout << "find(v1, 3)   = "; std::cout << tmp::deref<tmp::find<v1, tmp::integer<int, 3> >::type>::type::value << std::endl;
	std::cout << "find(v1, -32) = "; std::cout << tmp::deref<tmp::find<v1, tmp::integer<int, -32> >::type>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << std::boolalpha;
	std::cout << "Testing is_member" << std::endl;
	std::cout << "=================" << std::endl;
	std::cout << "v1          = "; print<v1>::eval();
	std::cout << "{3} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::integer<int, 3> >::type::value << std::endl;
	std::cout << "{2} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::integer<int, 2> >::type::value << std::endl;
	std::cout << "{0} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::integer<int, 0> >::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing serialise" << std::endl;
	std::cout << "=================" << std::endl;
	std::cout << "v1           = "; print<v1>::eval();
	std::cout << "v2           = "; print<v2>::eval();
	std::cout << "[v1, v2, v2] = "; print<tmp::serialise<tmp::vector<v1, v2, v2> >::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing unique" << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "v1         = "; print<v1>::eval();
	std::cout << "unique{v1} = "; print<tmp::unique<v1>::type>::eval();
	std::cout << std::endl;
	
	return 0;
}
