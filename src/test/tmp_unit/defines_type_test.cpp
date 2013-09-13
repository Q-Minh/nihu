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

#include "tmp/defines_type.hpp"
#include <iostream>

// metafunction declarations with macro
DEFINES_TYPE(default_version)
DEFINES_TYPE(special_version)

template <class C>
struct test
{
	typedef void default_version;
};

template <>
struct test<bool>
{
	typedef bool special_version;
};

int main(void)
{
	std::cout << std::boolalpha;

	std::cout << "Testing DEFAULT_VERSION" << std::endl;
	std::cout << "=======================" << std::endl;
	std::cout << "defines_default_version : " << defines_default_version<test<void> >::value << std::endl;
	std::cout << "defines_special_version : " << defines_special_version<test<void> >::value << std::endl;
	std::cout << std::endl;

	std::cout << "Testing SPECIAL_VERSION" << std::endl;
	std::cout << "=======================" << std::endl;
	std::cout << "defines_default_version : "  << defines_default_version<test<bool> >::value << std::endl;
	std::cout << "defines_special_version : " << defines_special_version<test<bool> >::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}
