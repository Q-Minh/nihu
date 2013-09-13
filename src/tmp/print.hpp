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

#ifndef PRINT_HPP_INCLUDED
#define PRINT_HPP_INCLUDED

#include <iostream>

template <class Seq>
struct print
{
	static void eval(void)
	{
		std::cout << tmp::deref<typename tmp::begin<Seq>::type>::type::value << ", ";
		print<typename tmp::pop_front<Seq>::type>::eval();
	}
};

template <>
struct print<tmp::vector<> >
{
	static void eval(void)
	{
		std::cout << std::endl;
	}
};

#endif // PRINT_HPP_INCLUDED
