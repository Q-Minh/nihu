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

#include "tmp/integer.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::int_<-31> _t;
	typedef tmp::int_<-4> _a;
	typedef tmp::int_< 5> _b;
	typedef tmp::int_< 8> _c;
	typedef tmp::int_<-9> _d;
	
	std::cout << std::boolalpha;

	std::cout << "Testing members with N = " << _t::value << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "N::value = " << _t::value << std::endl;
	std::cout << "N::next  = " << _t::next << std::endl;
	std::cout << "N::prev  = " << _t::prev << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing primitive metafunctions" << std::endl;
	std::cout << "===============================" << std::endl;
	std::cout << "next<N> = " << tmp::next<_t>::value << std::endl;
	std::cout << "prev<N> = " << tmp::prev<_t>::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing binary operations" << std::endl;
	std::cout << "=========================" << std::endl;
	std::cout << "a = " << _a::value << ", b = " << _b::value << std::endl;
	std::cout << "a + b = " << tmp::plus<_a, _b>::value << std::endl;
	std::cout << "a - b = " << tmp::minus<_a, _b>::value << std::endl;
	std::cout << "b - a = " << tmp::minus<_b, _a>::value << std::endl; 
	std::cout << "a * b = " << tmp::mul<_a, _b>::value << std::endl;
	std::cout << "a > b : " << tmp::greater<_a, _b>::value << std::endl;
	std::cout << "a < b : " << tmp::less<_a, _b>::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing extrema" << std::endl;
	std::cout << "===============" << std::endl;
	std::cout << "a = " << _a::value << ", b = " << _b::value << ", c = " << _c::value << ", d = " << _d::value << std::endl;
	std::cout << "max(a, b, c, d) = " << tmp::max_<_a, _b, _c, _d>::value << std::endl;
	std::cout << "min(a, b, c, d) = " << tmp::min_<_a, _b, _c, _d>::value << std::endl;
	
	return 0;
}
