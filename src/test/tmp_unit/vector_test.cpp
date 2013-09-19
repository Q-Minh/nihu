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
#include "tmp/vector.hpp"
#include "tmp/print.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::vector<> v0;
	typedef tmp::vector<tmp::integer<int, 1>, tmp::integer<int, -3>, tmp::integer<int, 0> > v1;
	typedef tmp::vector<tmp::integer<int, -4>, tmp::integer<int, 26>, tmp::integer<int, 7>, tmp::integer<int, 13> > v2;
	
	std::cout << "Test vectors" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "v0 : "; print<v0>::eval();
	std::cout << "v1 : "; print<v1>::eval();
	std::cout << "v2 : "; print<v2>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing at" << std::endl;
	std::cout << "==========" << std::endl;
	std::cout << "at<v1, 1> = " << tmp::at<v1, tmp::integer<int, 1> >::type::value << std::endl;
	std::cout << "at<v2, 3> = " << tmp::at<v2, tmp::integer<int, 3> >::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing pop and push" << std::endl;
	std::cout << "====================" << std::endl;
	std::cout << "push_back(v0, 16)   = "; print<tmp::push_back<v0, tmp::integer<int, 16> >::type>::eval();
	std::cout << "push_front(v2, -12) = "; print<tmp::push_front<v2, tmp::integer<int, -12> >::type>::eval();
	std::cout << "pop_front(v1)       = "; print<tmp::pop_front<v1>::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing size" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "size(v0) = " << tmp::size<v0>::value << std::endl;
	std::cout << std::endl;
	
	// Test vector iterators
	std::cout << "Testing vector iterators" << std::endl;
	std::cout << "========================" << std::endl;
	std::cout << "*begin(v1)       = " << tmp::deref<tmp::begin<v1>::type>::type::value << std::endl;
	std::cout << "*next(begin(v1)) = " << tmp::deref<tmp::next<tmp::begin<v1>::type>::type>::type::value << std::endl;
	std::cout << "*prev(end(v1))   = " << tmp::deref<tmp::prev<tmp::end<v1>::type>::type>::type::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}
