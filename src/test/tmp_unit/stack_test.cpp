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

#include "tmp/stack.hpp"
#include "tmp/algorithm.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::empty_stack myStack_0;
	typedef tmp::push_front<myStack_0, tmp::int_<1> >::type myStack_1;
	typedef tmp::push_front<myStack_1, tmp::int_<2> >::type myStack_2;
	typedef tmp::pop_front<myStack_2>::type myStack_3;
	
	std::cout << "Testing stack and operations" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "size(<>)                  = " << tmp::size<myStack_0>::type::value << std::endl;
	std::cout << "size(push(<>))            = " << tmp::size<myStack_1>::type::value << std::endl;
	std::cout << "size(push(push(<>)))      = " << tmp::size<myStack_2>::type::value << std::endl;
	std::cout << "size(pop(push(push(<>)))) = " << tmp::size<myStack_3>::type::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}

