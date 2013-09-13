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

#include "tmp/bool.hpp"

#include <type_traits>
#include <iostream>

int main(void)
{
	typedef std::true_type _t;
	typedef std::false_type _f;
	
	std::cout << std::boolalpha;
	
	std::cout << "Testing negation" << std::endl;
	std::cout << "================" << std::endl;
	std::cout << "!true  = " << tmp::not_<_t>::type::value << std::endl;
	std::cout << "!false = " << tmp::not_<_f>::type::value << std::endl; 
	std::cout << std::endl;
	
	std::cout << "Testing disjunction" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "true  | true  = " << tmp::or_<_t, _t>::type::value << std::endl;
	std::cout << "true  | false = " << tmp::or_<_t, _f>::type::value << std::endl;
	std::cout << "false | true  = " << tmp::or_<_f, _t>::type::value << std::endl;
	std::cout << "false | false = " << tmp::or_<_f, _f>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing conjunction" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "true  & true  = " << tmp::and_<_t, _t>::type::value << std::endl;
	std::cout << "true  & false = " << tmp::and_<_t, _f>::type::value << std::endl;
	std::cout << "false & true  = " << tmp::and_<_f, _t>::type::value << std::endl;
	std::cout << "false & false = " << tmp::and_<_f, _f>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << std::noboolalpha;
	
	std::cout << "Testing IF" << std::endl;
	std::cout << "==========" << std::endl;
	std::cout << "if_(true , 0, 1) = " << tmp::if_<_t, _f, _t >::type::value << std::endl;
	std::cout << "if_(false, 0, 1) = " << tmp::if_<_f, _f, _t >::type::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}
