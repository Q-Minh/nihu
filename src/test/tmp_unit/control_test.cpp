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

#include "tmp/algorithm.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

#include <type_traits>
#include <iostream>

template <class T>
struct call_tester
{
	struct type
	{
		void operator () (int const &, int const &)
		{
			std::cout << T::value << ' ';
		}
	};
};

template <class T, class U>
struct d_call_tester
{
	struct type
	{
		void operator () (int const &, int const &)
		{
			std::cout << "(" << T::value << ", " << U::value << ") ";
		}
	};
};

template <class T>
struct until_tester
{
	struct type
	{
		bool operator () (int const &)
		{
			std::cout << T::value << ' ';
			return T::value == 0;
		}
	};
};


int main(void)
{
	typedef tmp::vector<
		std::integral_constant<int, 1>,
		std::integral_constant<int, 2>
	> v1;
	
	typedef tmp::vector<
		std::integral_constant<int, -3>,
		std::integral_constant<int, 6>,
		std::integral_constant<int, 0>,
		std::integral_constant<int, -5>
	> v2;
	

	typedef call_tester<tmp::_1> call_test_t;
	typedef d_call_tester<tmp::_1, tmp::_2> d_call_test_t;
	typedef until_tester<tmp::_1> until_test_t;
	
	// Print test vectors
	std::cout << "Test vectors" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "v1 = [1, 2]" << std::endl;
	std::cout << "v2 = [-3 6 0 -5]" << std::endl;
	std::cout << std::endl;
	
	// Test call_each
	std::cout << "Testing call_each" << std::endl;
	std::cout << "=================" << std::endl;
	
	std::cout << "v1 = ";
	tmp::call_each<v1, call_test_t, int const &, int const &>(0, 0);
	std::cout << std::endl;
	
	std::cout << "v2 = ";
	tmp::call_each<v2, call_test_t, int const &, int const &>(0, 0);
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Test d_call_each
	std::cout << "Testing d_call_each" << std::endl;
	std::cout << "===================" << std::endl;
	
	std::cout << "v1 x v2 = ";
	tmp::d_call_each<v1, v2, d_call_test_t, int const &, int const&>(0, 0);
	std::cout << std::endl;
	
	std::cout << "v2 x v1 = ";
	tmp::d_call_each<v2, v1, d_call_test_t, int const &, int const&>(0, 0);
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Test call_until
	std::cout << std::boolalpha;
	std::cout << "Testing call_until" << std::endl;
	std::cout << "==================" << std::endl;
	std::cout << "v1 == 0? : ";
	std::cout << tmp::call_until<v1, until_test_t, int const&>(0) << std::endl;
	std::cout << "v2 == 0? : ";
	std::cout << tmp::call_until<v2, until_test_t, int const&>(0) << std::endl;
	std::cout << std::endl; 
	
	return 0;
}
